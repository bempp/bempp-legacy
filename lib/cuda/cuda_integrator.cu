// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "cuda_integrator.hpp"
#include "cuda_grid.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/shapeset.hpp"
#include "../fiber/basis_data.hpp"
#include "../fiber/scalar_traits.hpp"

#include <chrono>
#include <thrust/device_vector.h>

namespace Fiber {

template <typename BasisFunctionType, typename KernelType, typename ResultType>
CudaIntegrator<BasisFunctionType, KernelType, ResultType>::CudaIntegrator(
    const Matrix<CoordinateType> &localTestQuadPoints,
    const Matrix<CoordinateType> &localTrialQuadPoints,
    const std::vector<CoordinateType> &testQuadWeights,
    const std::vector<CoordinateType> &trialQuadWeights,
    const Shapeset<BasisFunctionType> &testShapeset,
    const Shapeset<BasisFunctionType> &trialShapeset,
    shared_ptr<Bempp::CudaGrid> testGrid,
    shared_ptr<Bempp::CudaGrid> trialGrid,
    const CollectionOfKernels<KernelType> &kernels,
    bool cacheElemData)
    : m_localTestQuadPoints(localTestQuadPoints),
      m_localTrialQuadPoints(localTrialQuadPoints),
      m_testGrid(testGrid), m_trialGrid(trialGrid),
      m_kernels(kernels), m_cacheElemData(cacheElemData) {

  const unsigned int testPointCount = localTestQuadPoints.cols();
  const unsigned int trialPointCount = localTrialQuadPoints.cols();

  if (testPointCount != testQuadWeights.size())
    throw std::invalid_argument(
        "CudaIntegrator::CudaIntegrator(): "
        "numbers of test points and weights do not match");
  if (trialPointCount != trialQuadWeights.size())
    throw std::invalid_argument(
        "CudaIntegrator::CudaIntegrator(): "
        "numbers of trial points and weights do not match");

  // Copy numerical quadrature weights to device memory
  m_testQuadData.weights = thrust::device_malloc<CoordinateType>(testPointCount);
  thrust::copy(testQuadWeights.begin(), testQuadWeights.end(),
      m_testQuadData.weights);
  if (testQuadWeights == trialQuadWeights) {
    m_trialQuadData.weights = m_testQuadData.weights;
  } else {
    m_trialQuadData.weights = thrust::device_malloc<CoordinateType>(trialPointCount);
      thrust::copy(trialQuadWeights.begin(), trialQuadWeights.end(),
          m_trialQuadData.weights);
  }
  m_testQuadData.pointCount = testPointCount;
  m_trialQuadData.pointCount = trialPointCount;

  cudaMemcpyToSymbol(constTestQuadWeights, &testQuadWeights[0],
      testPointCount*sizeof(CoordinateType));
  cudaMemcpyToSymbol(constTrialQuadWeights, &trialQuadWeights[0],
      trialPointCount*sizeof(CoordinateType));

  // Evaluate shapesets and copy basis data to device memory
  m_testBasisData.dofCount = testShapeset.size();
  m_trialBasisData.dofCount = trialShapeset.size();
  BasisData<BasisFunctionType> testBasisData, trialBasisData;
  size_t testBasisDeps = 0, trialBasisDeps = 0;
  testBasisDeps |= VALUES;
  trialBasisDeps |= VALUES;
  testShapeset.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS,
                        testBasisData);
  trialShapeset.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS,
                         trialBasisData);
  m_testBasisData.values = thrust::device_malloc<CudaBasisFunctionType>(
      testBasisData.componentCount() *
      testBasisData.functionCount() *
      testBasisData.pointCount());
  m_trialBasisData.values = thrust::device_malloc<CudaBasisFunctionType>(
      trialBasisData.componentCount() *
      trialBasisData.functionCount() *
      trialBasisData.pointCount());
  thrust::copy(testBasisData.values.begin(), testBasisData.values.end(),
      m_testBasisData.values);
  thrust::copy(trialBasisData.values.begin(), trialBasisData.values.end(),
      m_trialBasisData.values);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
CudaIntegrator<BasisFunctionType, KernelType, ResultType>::~CudaIntegrator() {

  thrust::device_free(m_testQuadData.weights);
  if (m_testQuadData.weights.get() != m_testQuadData.weights.get())
    thrust::device_free(m_trialQuadData.weights);

  thrust::device_free(m_testBasisData.values);
  thrust::device_free(m_trialBasisData.values);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
void CudaIntegrator<BasisFunctionType, KernelType, ResultType>::integrate(
    std::vector<int>::iterator startElementPairTestIndices,
    std::vector<int>::iterator endElementPairTestIndices,
    std::vector<int>::iterator startElementPairTrialIndices,
    std::vector<int>::iterator endElementPairTrialIndices,
    typename thrust::host_vector<ResultType>::iterator startResult) {

  const int testDofCount = m_testBasisData.dofCount;
  const int trialDofCount = m_trialBasisData.dofCount;

  const unsigned int testPointCount = m_localTestQuadPoints.cols();
  const unsigned int trialPointCount = m_localTrialQuadPoints.cols();

  const unsigned int testPointDim = m_localTestQuadPoints.rows();
  const unsigned int trialPointDim = m_localTrialQuadPoints.rows();

  if (testPointDim != 2 || trialPointDim != 2)
    throw std::invalid_argument("CudaIntegrator::integrate(): "
                                "only valid for two-dimensional local points");

  const unsigned int geometryPairCount =
      endElementPairTestIndices-startElementPairTestIndices;

  if (endElementPairTestIndices-startElementPairTestIndices
      != endElementPairTrialIndices-startElementPairTrialIndices)
    throw std::invalid_argument(
        "CudaIntegrator::integrate(): "
        "arrays 'elementPairTestIndices' and 'elementPairTrialIndices' must "
        "have the same number of elements");

  if (testPointCount == 0 || trialPointCount == 0 || geometryPairCount == 0)
    return;
  // TODO: in the (pathological) case that pointCount == 0 but
  // geometryPairCount != 0, set elements of result to 0.

  // Measure time of the GPU execution (CUDA event based)
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  // Copy element pair indices to device memory
  thrust::device_vector<int> d_elementPairTestIndices(
      startElementPairTestIndices, endElementPairTestIndices);
  thrust::device_vector<int> d_elementPairTrialIndices(
      startElementPairTrialIndices, endElementPairTrialIndices);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  float elapsedTimeIndexCopy;
  cudaEventElapsedTime(&elapsedTimeIndexCopy, start, stop);
  std::cout << "Time for thrust::copy(elementPairIndices, 'HtD') "
    << elapsedTimeIndexCopy << " ms" << std::endl;

  // Allocate device memory for the result
  thrust::device_vector<CudaResultType> d_result(
      geometryPairCount * testDofCount * trialDofCount);

  cudaEventRecord(start, 0);

  if (m_cacheElemData == true) {

    // Setup geometry data for selected elements on the device
    m_testGrid->setupGeometry();
    m_trialGrid->setupGeometry();

    ElemData<CoordinateType> testElemData, trialElemData;

    thrust::device_vector<CoordinateType> d_testGeomData, d_trialGeomData;
    thrust::device_vector<CoordinateType> d_testNormals, d_trialNormals;
    thrust::device_vector<CoordinateType> d_testIntegrationElements, d_trialIntegrationElements;

    // Calculate global points on the device
    m_testGrid->local2global(m_localTestQuadPoints, d_testGeomData);
    testElemData.geomData = d_testGeomData.data();
    if (m_testGrid.get() == m_trialGrid.get() &&
        m_localTestQuadPoints == m_localTrialQuadPoints) {
      trialElemData.geomData = testElemData.geomData;
    } else {
      m_trialGrid->local2global(m_localTrialQuadPoints, d_trialGeomData);
      trialElemData.geomData = d_trialGeomData.data();
    }

    m_testGrid->calculateNormalsAndIntegrationElements(
        d_testNormals, d_testIntegrationElements);
    testElemData.normals = d_testNormals.data();
    testElemData.integrationElements = d_testIntegrationElements.data();
    testElemData.activeElemCount = d_testIntegrationElements.size();
    if (m_testGrid.get() == m_trialGrid.get()) {
      trialElemData.normals = testElemData.normals;
      trialElemData.integrationElements = testElemData.integrationElements;
      trialElemData.activeElemCount = testElemData.activeElemCount;
    } else {
      m_trialGrid->calculateNormalsAndIntegrationElements(
          d_trialNormals, d_trialIntegrationElements);
      trialElemData.normals = d_trialNormals.data();
      trialElemData.integrationElements = d_trialIntegrationElements.data();
      trialElemData.activeElemCount = d_trialIntegrationElements.size();
    }

    // Each thread is working on one pair of elements
    thrust::counting_iterator<int> iter(0);
    thrust::for_each(iter, iter+geometryPairCount,
        CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>(
            d_elementPairTestIndices.data(),
            d_elementPairTrialIndices.data(),
            m_testQuadData, m_trialQuadData,
            m_testBasisData, m_trialBasisData,
            testElemData, trialElemData,
            d_result.data()
            ));
  } else {

    // Get raw geometry data
    RawGeometryData<double> testRawGeometryData, trialRawGeometryData;
    m_testGrid->getRawGeometryData(
        testRawGeometryData.vtxCount, testRawGeometryData.elemCount,
        testRawGeometryData.vertices, testRawGeometryData.elementCorners);
    m_trialGrid->getRawGeometryData(
        trialRawGeometryData.vtxCount, trialRawGeometryData.elemCount,
        trialRawGeometryData.vertices, trialRawGeometryData.elementCorners);

    thrust::host_vector<CoordinateType> h_testFun0(testPointCount);
    thrust::host_vector<CoordinateType> h_testFun1(testPointCount);
    thrust::host_vector<CoordinateType> h_testFun2(testPointCount);
    thrust::host_vector<CoordinateType> h_trialFun0(trialPointCount);
    thrust::host_vector<CoordinateType> h_trialFun1(trialPointCount);
    thrust::host_vector<CoordinateType> h_trialFun2(trialPointCount);

    // Evaluate geometrical shape function values
    for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
      const CoordinateType r = m_localTestQuadPoints(0,testPoint);
      const CoordinateType s = m_localTestQuadPoints(1,testPoint);
      h_testFun0[testPoint] = 1.0 - r - s;
      h_testFun1[testPoint] = r;
      h_testFun2[testPoint] = s;
    }
    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      const CoordinateType r = m_localTrialQuadPoints(0,trialPoint);
      const CoordinateType s = m_localTrialQuadPoints(1,trialPoint);
      h_trialFun0[trialPoint] = 1.0 - r - s;
      h_trialFun1[trialPoint] = r;
      h_trialFun2[trialPoint] = s;
    }

    // Copy data to device
    thrust::device_vector<CoordinateType> d_testFun0 = h_testFun0;
    thrust::device_vector<CoordinateType> d_testFun1 = h_testFun1;
    thrust::device_vector<CoordinateType> d_testFun2 = h_testFun2;
    thrust::device_vector<CoordinateType> d_trialFun0 = h_trialFun0;
    thrust::device_vector<CoordinateType> d_trialFun1 = h_trialFun1;
    thrust::device_vector<CoordinateType> d_trialFun2 = h_trialFun2;

    // Collect geometrical shape function data
    GeomShapeFunData<CoordinateType> testGeomShapeFunData, trialGeomShapeFunData;
    testGeomShapeFunData.fun0 = d_testFun0.data();
    testGeomShapeFunData.fun1 = d_testFun1.data();
    testGeomShapeFunData.fun2 = d_testFun2.data();
    trialGeomShapeFunData.fun0 = d_trialFun0.data();
    trialGeomShapeFunData.fun1 = d_trialFun1.data();
    trialGeomShapeFunData.fun2 = d_trialFun2.data();

    // Copy geometrical shape function data to constant device memory
    cudaMemcpyToSymbol(constTestGeomShapeFun0, h_testFun0.data(),
        testPointCount*sizeof(CoordinateType));
    cudaMemcpyToSymbol(constTestGeomShapeFun1, h_testFun1.data(),
        testPointCount*sizeof(CoordinateType));
    cudaMemcpyToSymbol(constTestGeomShapeFun2, h_testFun2.data(),
        testPointCount*sizeof(CoordinateType));
    cudaMemcpyToSymbol(constTrialGeomShapeFun0, h_trialFun0.data(),
        trialPointCount*sizeof(CoordinateType));
    cudaMemcpyToSymbol(constTrialGeomShapeFun1, h_trialFun1.data(),
        trialPointCount*sizeof(CoordinateType));
    cudaMemcpyToSymbol(constTrialGeomShapeFun2, h_trialFun2.data(),
        trialPointCount*sizeof(CoordinateType));

    // Each thread is working on one pair of elements
    thrust::counting_iterator<int> iter(0);
    thrust::for_each(iter, iter+geometryPairCount,
        CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorNonCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>(
            d_elementPairTestIndices.data(),
            d_elementPairTrialIndices.data(),
            m_testQuadData, m_trialQuadData,
            m_testBasisData, m_trialBasisData,
            testRawGeometryData, trialRawGeometryData,
            testGeomShapeFunData, trialGeomShapeFunData,
            d_result.data()
            ));

//    unsigned int blockSze = 512;
//    dim3 blockSize(blockSze,1,1);
//    unsigned int gridSze = std::max(
//        static_cast<unsigned int>((geometryPairCount-1)/blockSze+1),
//        static_cast<unsigned int>(1));
//    dim3 gridSize(gridSze,1,1);
//    std::cout << "gridSize = " << gridSize.x << std::endl;
//
//    RawCudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorNonCached<
//        CudaBasisFunctionType, CudaKernelType, CudaResultType>
//        <<<gridSize,blockSize>>>(
//        geometryPairCount,
//        thrust::raw_pointer_cast(d_elementPairTestIndices.data()),
//        thrust::raw_pointer_cast(d_elementPairTrialIndices.data()),
//        m_testQuadData.pointCount, m_trialQuadData.pointCount,
//        m_testBasisData.dofCount, thrust::raw_pointer_cast(m_testBasisData.values),
//        m_trialBasisData.dofCount, thrust::raw_pointer_cast(m_trialBasisData.values),
//        testRawGeometryData.elemCount, testRawGeometryData.vtxCount,
//        thrust::raw_pointer_cast(testRawGeometryData.vertices),
//        thrust::raw_pointer_cast(testRawGeometryData.elementCorners),
//        trialRawGeometryData.elemCount, trialRawGeometryData.vtxCount,
//        thrust::raw_pointer_cast(trialRawGeometryData.vertices),
//        thrust::raw_pointer_cast(trialRawGeometryData.elementCorners),
//        thrust::raw_pointer_cast(d_result.data()));
  }

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  float elapsedTimeIntegral;
  cudaEventElapsedTime(&elapsedTimeIntegral , start, stop);
  std::cout << "Time for CudaEvaluateIntegral() is "
    << elapsedTimeIntegral << " ms" << std::endl;

  cudaEventRecord(start, 0);

  // Copy result back to host memory
  thrust::copy(d_result.begin(), d_result.end(), startResult);

//  std::cout << "h_result = " << std::endl;
//  for (int i = 0; i < h_result.size(); ++i) {
//    std::cout << h_result[i] << " " << std::flush;
//  }
//  std::cout << std::endl;

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  float elapsedTimeResultCopy;
  cudaEventElapsedTime(&elapsedTimeResultCopy, start, stop);
  std::cout << "Time for thrust::copy(result, 'DtH') "
    << elapsedTimeResultCopy << " ms" << std::endl;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(CudaIntegrator);

} // namespace Fiber
