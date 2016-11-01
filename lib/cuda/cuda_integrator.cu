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

#include <complex>
#include <chrono>
#include <thrust/device_vector.h>

namespace Fiber {

template <typename BasisFunctionType, typename ResultType>
CudaIntegrator<BasisFunctionType, ResultType>::CudaIntegrator(
    const Matrix<CoordinateType> &localTestQuadPoints,
    const Matrix<CoordinateType> &localTrialQuadPoints,
    const std::vector<CoordinateType> &testQuadWeights,
    const std::vector<CoordinateType> &trialQuadWeights,
    const Shapeset<BasisFunctionType> &testShapeset,
    const Shapeset<BasisFunctionType> &trialShapeset,
    shared_ptr<Bempp::CudaGrid> testGrid,
    shared_ptr<Bempp::CudaGrid> trialGrid,
    bool cacheElemData)
    : m_localTestQuadPoints(localTestQuadPoints),
      m_localTrialQuadPoints(localTrialQuadPoints),
      m_testGrid(testGrid), m_trialGrid(trialGrid),
      m_cacheElemData(cacheElemData) {

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
  testShapeset.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS,
                        testBasisData);
  trialShapeset.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS,
                         trialBasisData);
  m_testBasisData.basisData = thrust::device_malloc<BasisFunctionType>(
      testBasisData.componentCount() *
      testBasisData.functionCount() *
      testBasisData.pointCount());
  m_trialBasisData.basisData = thrust::device_malloc<BasisFunctionType>(
      trialBasisData.componentCount() *
      trialBasisData.functionCount() *
      trialBasisData.pointCount());
  thrust::copy(testBasisData.values.begin(), testBasisData.values.end(),
      m_testBasisData.basisData);
  thrust::copy(trialBasisData.values.begin(), trialBasisData.values.end(),
      m_trialBasisData.basisData);
}

template <typename BasisFunctionType, typename ResultType>
CudaIntegrator<BasisFunctionType, ResultType>::~CudaIntegrator() {

  thrust::device_free(m_testQuadData.weights);
  if (m_testQuadData.weights.get() != m_testQuadData.weights.get())
    thrust::device_free(m_trialQuadData.weights);

  thrust::device_free(m_testBasisData.basisData);
  thrust::device_free(m_trialBasisData.basisData);
}

template <typename BasisFunctionType, typename ResultType>
void CudaIntegrator<BasisFunctionType, ResultType>::integrate(
    const std::vector<int> &elementPairTestIndices,
    const std::vector<int> &elementPairTrialIndices,
    std::vector<Matrix<ResultType>*> &result) {

  std::cout << "Hello, this is CudaIntegrator::integrate()!" << std::endl;

  const unsigned int testPointCount = m_localTestQuadPoints.cols();
  const unsigned int trialPointCount = m_localTrialQuadPoints.cols();

  const unsigned int testPointDim = m_localTestQuadPoints.rows();
  const unsigned int trialPointDim = m_localTrialQuadPoints.rows();

  if (testPointDim != 2 || trialPointDim != 2)
    throw std::invalid_argument("CudaIntegrator::integrate(): "
                                "only valid for two-dimensional local points");

  const unsigned int geometryPairCount = elementPairTestIndices.size();

  if (elementPairTestIndices.size() != elementPairTrialIndices.size())
    throw std::invalid_argument(
        "CudaIntegrator::integrate(): "
        "arrays 'elementPairTestIndices' and 'elementPairTrialIndices' must "
        "have the same number of elements");

  if (result.size() != geometryPairCount)
    throw std::invalid_argument(
        "CudaIntegrator::integrate(): "
        "arrays 'result' and 'elementPairIndices' must have the same number "
        "of elements");

  if (testPointCount == 0 || trialPointCount == 0 || geometryPairCount == 0)
    return;
  // TODO: in the (pathological) case that pointCount == 0 but
  // geometryPairCount != 0, set elements of result to 0.

  const int testDofCount = m_testBasisData.dofCount;
  const int trialDofCount = m_trialBasisData.dofCount;

  for (size_t i = 0; i < result.size(); ++i) {
    assert(result[i]);
    result[i]->resize(testDofCount, trialDofCount);
  }

  // Copy element pair indices to device memory
  thrust::device_vector<int> d_elementPairTestIndices(elementPairTestIndices);
  thrust::device_vector<int> d_elementPairTrialIndices(elementPairTrialIndices);

  // Allocate device memory for the result
  thrust::device_vector<CoordinateType> d_result(
      geometryPairCount * testDofCount * trialDofCount);

  // Measure time of the GPU execution (CUDA event based)
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
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
    typedef CoordinateType KernelType;
    thrust::counting_iterator<int> iter(0);
    thrust::for_each(iter, iter+geometryPairCount,
        EvaluateLaplace3dSingleLayerPotentialIntegralFunctorCached<
        BasisFunctionType, KernelType, CoordinateType>(
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
    typedef CoordinateType KernelType;
    thrust::counting_iterator<int> iter(0);
    thrust::for_each(iter, iter+geometryPairCount,
        EvaluateLaplace3dDoubleLayerPotentialIntegralFunctorNonCached<
        BasisFunctionType, KernelType, CoordinateType>(
            d_elementPairTestIndices.data(),
            d_elementPairTrialIndices.data(),
            m_testQuadData, m_trialQuadData,
            m_testBasisData, m_trialBasisData,
            testRawGeometryData, trialRawGeometryData,
            testGeomShapeFunData, trialGeomShapeFunData,
            d_result.data()
            ));
  }

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  float elapsedTimeIntegral;
  cudaEventElapsedTime(&elapsedTimeIntegral , start, stop);
  std::cout << "Time for integral evaluation is "
    << elapsedTimeIntegral << " ms" << std::endl;

  cudaEventRecord(start, 0);

  // Copy result back to host memory
  thrust::host_vector<CoordinateType> h_result = d_result;

//  std::cout << "h_result = " << std::endl;
//  for (int i = 0; i < h_result.size(); ++i) {
//    std::cout << h_result[i] << " " << std::flush;
//  }
//  std::cout << std::endl;

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  float elapsedTimeResultCopy;
  cudaEventElapsedTime(&elapsedTimeResultCopy, start, stop);
  std::cout << "Time for result copy is "
    << elapsedTimeResultCopy << " ms" << std::endl;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  // Assemble result
  for (int geometryPair = 0; geometryPair < geometryPairCount; ++geometryPair) {
    for (int testDof = 0; testDof < testDofCount; ++testDof) {
      for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
        (*result[geometryPair])(testDof, trialDof) =
            h_result[geometryPair * testDofCount * trialDofCount
                   + testDof * trialDofCount
                   + trialDof];
      }
    }
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time for local result assembly = "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << " ms" << std::endl;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(CudaIntegrator);

} // namespace Fiber
