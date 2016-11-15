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

namespace Fiber {

template <typename BasisFunctionType, typename KernelType, typename ResultType>
CudaIntegrator<BasisFunctionType, KernelType, ResultType>::CudaIntegrator(
    const Matrix<CoordinateType> &localTestQuadPoints,
    const Matrix<CoordinateType> &localTrialQuadPoints,
    const std::vector<CoordinateType> &testQuadWeights,
    const std::vector<CoordinateType> &trialQuadWeights,
    const Shapeset<BasisFunctionType> &testShapeset,
    const Shapeset<BasisFunctionType> &trialShapeset,
    shared_ptr<Bempp::CudaGrid<CoordinateType>> testGrid,
    shared_ptr<Bempp::CudaGrid<CoordinateType>> trialGrid,
    const std::vector<int> &testIndices, const std::vector<int> &trialIndices,
    const CollectionOfKernels<KernelType> &kernels,
    const int deviceId, const Bempp::CudaOptions &cudaOptions)
    : m_localTestQuadPoints(localTestQuadPoints),
      m_localTrialQuadPoints(localTrialQuadPoints),
      m_testGrid(testGrid), m_trialGrid(trialGrid),
      m_kernels(kernels), m_deviceId(deviceId), m_cudaOptions(cudaOptions) {

  const unsigned int testIndexCount = testIndices.size();
  const unsigned int trialIndexCount = trialIndices.size();

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

  cudaSetDevice(m_deviceId);

  // Copy element indices to device memory
  m_d_testIndices.resize(testIndexCount);
  m_d_trialIndices.resize(trialIndexCount);
  m_d_testIndices.assign(testIndices.begin(), testIndices.end());
  m_d_trialIndices.assign(trialIndices.begin(), trialIndices.end());

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

  if (cudaOptions.isElementDataCachingEnabled() == true) {

    // Setup geometry data for selected elements on the device
    m_testGrid->setupGeometry();
    m_trialGrid->setupGeometry();

    // Calculate global points on the device
    m_testGrid->local2global(m_localTestQuadPoints, m_d_testGeomData);
    m_testElemData.geomData = m_d_testGeomData.data();
    if (m_testGrid.get() == m_trialGrid.get() &&
        m_localTestQuadPoints == m_localTrialQuadPoints) {
      m_trialElemData.geomData = m_testElemData.geomData;
    } else {
      m_trialGrid->local2global(m_localTrialQuadPoints, m_d_trialGeomData);
      m_trialElemData.geomData = m_d_trialGeomData.data();
    }

    m_testGrid->calculateNormalsAndIntegrationElements(
        m_d_testNormals, m_d_testIntegrationElements);
    m_testElemData.normals = m_d_testNormals.data();
    m_testElemData.integrationElements = m_d_testIntegrationElements.data();
    m_testElemData.activeElemCount = m_d_testIntegrationElements.size();
    if (m_testGrid.get() == m_trialGrid.get()) {
      m_trialElemData.normals = m_testElemData.normals;
      m_trialElemData.integrationElements = m_testElemData.integrationElements;
      m_trialElemData.activeElemCount = m_testElemData.activeElemCount;
    } else {
      m_trialGrid->calculateNormalsAndIntegrationElements(
          m_d_trialNormals, m_d_trialIntegrationElements);
      m_trialElemData.normals = m_d_trialNormals.data();
      m_trialElemData.integrationElements = m_d_trialIntegrationElements.data();
      m_trialElemData.activeElemCount = m_d_trialIntegrationElements.size();
    }

  } else {

    // Get raw geometry data
    m_testGrid->getRawGeometryData(
        m_testRawGeometryData.vtxCount, m_testRawGeometryData.elemCount,
        m_testRawGeometryData.vertices, m_testRawGeometryData.elementCorners);
    m_trialGrid->getRawGeometryData(
        m_trialRawGeometryData.vtxCount, m_trialRawGeometryData.elemCount,
        m_trialRawGeometryData.vertices, m_trialRawGeometryData.elementCorners);

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
  }

  // Create CUDA streams
  const int streamCount = cudaOptions.streamCount();
  m_streams.resize(streamCount);
  for (int stream = 0; stream < streamCount; ++stream)
    cudaStreamCreate(&m_streams[stream]);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
CudaIntegrator<BasisFunctionType, KernelType, ResultType>::~CudaIntegrator() {

  cudaSetDevice(m_deviceId);

  thrust::device_free(m_testQuadData.weights);
  if (m_testQuadData.weights.get() != m_testQuadData.weights.get())
    thrust::device_free(m_trialQuadData.weights);

  thrust::device_free(m_testBasisData.values);
  thrust::device_free(m_trialBasisData.values);

  for (int stream = 0; stream < m_streams.size(); ++stream)
    cudaStreamDestroy(m_streams[stream]);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
void CudaIntegrator<BasisFunctionType, KernelType, ResultType>::integrate(
    const int elemPairIndexBegin, const int elemPairIndexEnd,
    ResultType *result) {

  const int testDofCount = m_testBasisData.dofCount;
  const int trialDofCount = m_trialBasisData.dofCount;

  const unsigned int testPointCount = m_localTestQuadPoints.cols();
  const unsigned int trialPointCount = m_localTrialQuadPoints.cols();

  const unsigned int testPointDim = m_localTestQuadPoints.rows();
  const unsigned int trialPointDim = m_localTrialQuadPoints.rows();

  if (testPointDim != 2 || trialPointDim != 2)
    throw std::invalid_argument("CudaIntegrator::integrate(): "
                                "only valid for two-dimensional local points");

  const unsigned int elemPairCount = elemPairIndexEnd - elemPairIndexBegin;

  const unsigned int streamCount = m_streams.size();
  const unsigned int streamSize =
      static_cast<unsigned int>((elemPairCount-1) / streamCount + 1);

  if (testPointCount == 0 || trialPointCount == 0 || elemPairCount == 0)
    return;
  // TODO: in the (pathological) case that pointCount == 0 but
  // geometryPairCount != 0, set elements of result to 0.

  cudaSetDevice(m_deviceId);

  // Allocate device memory for the result
  thrust::device_vector<CudaResultType> d_result(
      elemPairCount * testDofCount * trialDofCount);

  unsigned int blockSze = m_cudaOptions.blockSize();
  dim3 blockSize(blockSze,1,1);

  if (m_cudaOptions.isElementDataCachingEnabled() == true) {

    for (int stream = 0; stream < streamCount; ++stream) {

      int offset = stream * streamSize;

      int actualStreamSize;
      if (stream == streamCount-1) {
        actualStreamSize = elemPairCount - (streamCount-1) * streamSize;
      } else {
        actualStreamSize = streamSize;
      }

//      // Each thread is working on one pair of elements
//      thrust::counting_iterator<int> iter(0);
//      thrust::for_each(iter, iter+geometryPairCount,
//          CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorCached<
//          CudaBasisFunctionType, CudaKernelType, CudaResultType>(
//              d_elementPairTestIndices.data(),
//              d_elementPairTrialIndices.data(),
//              m_testQuadData, m_trialQuadData,
//              m_testBasisData, m_trialBasisData,
//              testElemData, trialElemData,
//              d_result.data()
//              ));

      unsigned int gridSze =
          static_cast<unsigned int>((actualStreamSize-1)/blockSze+1);
      dim3 gridSize(gridSze,1,1);

      // Launch kernel
      RawCudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorCached<
          CudaBasisFunctionType, CudaKernelType, CudaResultType>
          <<<gridSize, blockSize, 0, m_streams[stream]>>>(
          elemPairIndexBegin, m_d_trialIndices.size(),
          actualStreamSize, offset,
          thrust::raw_pointer_cast(m_d_testIndices.data()),
          thrust::raw_pointer_cast(m_d_trialIndices.data()),
          m_testQuadData.pointCount, m_trialQuadData.pointCount,
          m_testBasisData.dofCount,
          thrust::raw_pointer_cast(m_testBasisData.values),
          m_trialBasisData.dofCount,
          thrust::raw_pointer_cast(m_trialBasisData.values),
          m_testElemData.activeElemCount,
          thrust::raw_pointer_cast(m_testElemData.geomData),
          thrust::raw_pointer_cast(m_testElemData.normals),
          thrust::raw_pointer_cast(m_testElemData.integrationElements),
          m_trialElemData.activeElemCount,
          thrust::raw_pointer_cast(m_trialElemData.geomData),
          thrust::raw_pointer_cast(m_trialElemData.normals),
          thrust::raw_pointer_cast(m_trialElemData.integrationElements),
          thrust::raw_pointer_cast(d_result.data()+offset*testDofCount*trialDofCount));

      // Copy result to host memory
      cudaMemcpyAsync(
          &result[offset*testDofCount*trialDofCount],
          thrust::raw_pointer_cast(d_result.data()+offset*testDofCount*trialDofCount),
          actualStreamSize*testDofCount*trialDofCount*sizeof(ResultType),
          cudaMemcpyDeviceToHost, m_streams[stream]);
    }

  } else {

    for (int stream = 0; stream < streamCount; ++stream) {

      int offset = stream * streamSize;

      int actualStreamSize;
      if (stream == streamCount-1) {
        actualStreamSize = elemPairCount - (streamCount-1) * streamSize;
      } else {
        actualStreamSize = streamSize;
      }

      // Each thread is working on one pair of elements
//      thrust::counting_iterator<int> iter(0);
//      thrust::for_each(iter, iter+geometryPairCount,
//          CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorNonCached<
//          CudaBasisFunctionType, CudaKernelType, CudaResultType>(
//              d_elementPairTestIndices.data(),
//              d_elementPairTrialIndices.data(),
//              m_testQuadData, m_trialQuadData,
//              m_testBasisData, m_trialBasisData,
//              testRawGeometryData, trialRawGeometryData,
//              testGeomShapeFunData, trialGeomShapeFunData,
//              d_result.data()
//              ));

      unsigned int gridSze =
          static_cast<unsigned int>((actualStreamSize-1)/blockSze+1);
      dim3 gridSize(gridSze,1,1);

      // Launch kernel
      RawCudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorNonCached<
          CudaBasisFunctionType, CudaKernelType, CudaResultType>
          <<<gridSize, blockSize, 0, m_streams[stream]>>>(
          elemPairIndexBegin, m_d_trialIndices.size(),
          actualStreamSize, offset,
          thrust::raw_pointer_cast(m_d_testIndices.data()),
          thrust::raw_pointer_cast(m_d_trialIndices.data()),
          m_testQuadData.pointCount, m_trialQuadData.pointCount,
          m_testBasisData.dofCount, thrust::raw_pointer_cast(m_testBasisData.values),
          m_trialBasisData.dofCount, thrust::raw_pointer_cast(m_trialBasisData.values),
          m_testRawGeometryData.elemCount, m_testRawGeometryData.vtxCount,
          thrust::raw_pointer_cast(m_testRawGeometryData.vertices),
          thrust::raw_pointer_cast(m_testRawGeometryData.elementCorners),
          m_trialRawGeometryData.elemCount, m_trialRawGeometryData.vtxCount,
          thrust::raw_pointer_cast(m_trialRawGeometryData.vertices),
          thrust::raw_pointer_cast(m_trialRawGeometryData.elementCorners),
          thrust::raw_pointer_cast(d_result.data()+offset*testDofCount*trialDofCount));

      // Copy result to host memory
      cudaMemcpyAsync(
          &result[offset*testDofCount*trialDofCount],
          thrust::raw_pointer_cast(d_result.data()+offset*testDofCount*trialDofCount),
          actualStreamSize*testDofCount*trialDofCount*sizeof(ResultType),
          cudaMemcpyDeviceToHost, m_streams[stream]);
    }
  }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(CudaIntegrator);

} // namespace Fiber
