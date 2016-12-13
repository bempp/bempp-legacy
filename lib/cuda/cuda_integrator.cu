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

#include "cuda_evaluate_laplace_3d_single_layer_potential_integral.cuh"
#include "cuda_evaluate_laplace_3d_double_layer_potential_integral.cuh"
#include "cuda_evaluate_laplace_3d_adjoint_double_layer_potential_integral.cuh"
#include "cuda_evaluate_modified_helmholtz_3d_single_layer_potential_integral.cuh"
#include "cuda_evaluate_modified_helmholtz_3d_double_layer_potential_integral.cuh"
#include "cuda_evaluate_modified_helmholtz_3d_adjoint_double_layer_potential_integral.cuh"
#include "cuda_evaluate_modified_helmholtz_3d_hypersingular_integral.cuh"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/shapeset.hpp"
#include "../fiber/basis_data.hpp"
#include "../fiber/collection_of_kernels.hpp"

// Helper
namespace {

template <typename CudaBasisFunctionType, typename BasisFunctionType>
struct SetupBasisDataHelper {
  static void reorder(
      std::vector<CudaBasisFunctionType> &basisValues,
      std::vector<CudaBasisFunctionType> &basisDerivatives,
      const Fiber::BasisData<BasisFunctionType> &basisData,
      const unsigned int dofCount, const unsigned int pointCount,
      const unsigned int compCount,
      const int dof, const int point, const int comp) {
    throw std::invalid_argument(
        "CudaIntegrator::setupBasisData(): "
        "invalid BasisFunctionType");
  }
};

template <typename CudaBasisFunctionType>
struct SetupBasisDataHelper<CudaBasisFunctionType, float> {
  static void reorder(
      std::vector<CudaBasisFunctionType> &basisValues,
      std::vector<CudaBasisFunctionType> &basisDerivatives,
      const Fiber::BasisData<float> &basisData,
      const unsigned int dofCount, const unsigned int pointCount,
      const unsigned int compCount,
      const int dof, const int point, const int comp) {

    const unsigned int localCoordCount = 2;

    basisValues[dof * pointCount * compCount
                + point * compCount
                + comp] =
        basisData.values(comp, dof, point);

    for (unsigned int coordIndex = 0; coordIndex < localCoordCount; ++coordIndex) {
      basisDerivatives[coordIndex * dofCount * pointCount * compCount
                       + dof * pointCount * compCount
                       + point * compCount
                       + comp] =
          basisData.derivatives(comp, coordIndex, dof, point);
    }
  }
};

template <typename CudaBasisFunctionType>
struct SetupBasisDataHelper<CudaBasisFunctionType, double> {
  static void reorder(
      std::vector<CudaBasisFunctionType> &basisValues,
      std::vector<CudaBasisFunctionType> &basisDerivatives,
      const Fiber::BasisData<double> &basisData,
      const unsigned int dofCount, const unsigned int pointCount,
      const unsigned int compCount,
      const int dof, const int point, const int comp) {

    const unsigned int localCoordCount = 2;

    basisValues[dof * pointCount * compCount
                + point * compCount
                + comp] =
        basisData.values(comp, dof, point);

    for (unsigned int coordIndex = 0; coordIndex < localCoordCount; ++coordIndex) {
      basisDerivatives[coordIndex * dofCount * pointCount * compCount
                       + dof * pointCount * compCount
                       + point * compCount
                       + comp] =
          basisData.derivatives(comp, coordIndex, dof, point);
    }
  }
};

template <typename CudaKernelType, typename KernelType>
struct GetWaveNumberHelper {
  static void getWaveNumber(
      CudaKernelType &waveNumberReal, CudaKernelType &waveNumberImag,
      const Fiber::shared_ptr<
          const Fiber::CollectionOfKernels<KernelType>> &kernel) {
    throw std::invalid_argument(
        "CudaIntegrator::getWaveNumber(): "
        "invalid KernelType");
  }
};

template <typename CudaKernelType>
struct GetWaveNumberHelper<CudaKernelType, float> {
  static void getWaveNumber(
      CudaKernelType &waveNumberReal, CudaKernelType &waveNumberImag,
      const Fiber::shared_ptr<
          const Fiber::CollectionOfKernels<float>> &kernel) {
    waveNumberReal = static_cast<CudaKernelType>(0);
    waveNumberImag = static_cast<CudaKernelType>(0);
  }
};

template <typename CudaKernelType>
struct GetWaveNumberHelper<CudaKernelType, double> {
  static void getWaveNumber(
      CudaKernelType &waveNumberReal, CudaKernelType &waveNumberImag,
      const Fiber::shared_ptr<
          const Fiber::CollectionOfKernels<double>> &kernel) {
    waveNumberReal = static_cast<CudaKernelType>(0);
    waveNumberImag = static_cast<CudaKernelType>(0);
  }
};

template <typename CudaKernelType>
struct GetWaveNumberHelper<CudaKernelType, std::complex<float>> {
  static void getWaveNumber(
      CudaKernelType &waveNumberReal, CudaKernelType &waveNumberImag,
      const Fiber::shared_ptr<
          const Fiber::CollectionOfKernels<std::complex<float>>> &kernel) {
    waveNumberReal = kernel->waveNumber().real();
    waveNumberImag = kernel->waveNumber().imag();
  }
};

template <typename CudaKernelType>
struct GetWaveNumberHelper<CudaKernelType, std::complex<double>> {
  static void getWaveNumber(
      CudaKernelType &waveNumberReal, CudaKernelType &waveNumberImag,
      const Fiber::shared_ptr<
          const Fiber::CollectionOfKernels<std::complex<double>>> &kernel) {
    waveNumberReal = kernel->waveNumber().real();
    waveNumberImag = kernel->waveNumber().imag();
  }
};

} // namespace

namespace Fiber {

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
CudaIntegrator<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    CudaIntegrator(
        const Matrix<CoordinateType> &localTestQuadPoints,
        const Matrix<CoordinateType> &localTrialQuadPoints,
        const std::vector<CoordinateType> &testQuadWeights,
        const std::vector<CoordinateType> &trialQuadWeights,
        const Shapeset<BasisFunctionType> &testShapeset,
        const Shapeset<BasisFunctionType> &trialShapeset,
        const shared_ptr<Bempp::CudaGrid<CudaCoordinateType>> &testGrid,
        const shared_ptr<Bempp::CudaGrid<CudaCoordinateType>> &trialGrid,
        const std::vector<int> &testIndices, const std::vector<int> &trialIndices,
        const size_t maxActiveElemPairCount,
        const shared_ptr<const CollectionOfKernels<KernelType>> &kernel,
        const int deviceId, const Bempp::CudaOptions &cudaOptions)
        : m_testGrid(testGrid), m_trialGrid(trialGrid), m_kernelName(kernel->name()),
          m_deviceId(deviceId), m_cudaOptions(cudaOptions),
          m_isElementDataCachingEnabled(cudaOptions.isElementDataCachingEnabled()),
          m_isKernelDataCachingEnabled(cudaOptions.isKernelDataCachingEnabled()) {

  // Get wave number
  GetWaveNumberHelper<CudaKernelType, KernelType>::getWaveNumber(
      m_waveNumberReal, m_waveNumberImag, kernel);

  const unsigned int trialPointCount = localTrialQuadPoints.cols();
  const unsigned int testPointCount = localTestQuadPoints.cols();

  if (trialPointCount != trialQuadWeights.size())
    throw std::invalid_argument(
        "CudaIntegrator::CudaIntegrator(): "
        "numbers of trial points and weights do not match");
  if (testPointCount != testQuadWeights.size())
    throw std::invalid_argument(
        "CudaIntegrator::CudaIntegrator(): "
        "numbers of test points and weights do not match");

  if (trialPointCount > 6)
    throw std::invalid_argument(
        "CudaIntegrator::CudaIntegrator(): "
        "number of trial points too high in terms of device memory");
  if (testPointCount > 6)
    throw std::invalid_argument(
        "CudaIntegrator::CudaIntegrator(): "
        "number of test points too high in terms of device memory");

  cu_verify( cudaSetDevice(m_deviceId) );

  // Copy element indices to device memory
  m_d_trialIndices.resize(trialIndices.size());
  m_d_testIndices.resize(testIndices.size());
  m_d_trialIndices.assign(trialIndices.begin(), trialIndices.end());
  m_d_testIndices.assign(testIndices.begin(), testIndices.end());

  // Copy numerical quadrature weights to constant device memory
  m_trialQuadData.pointCount = trialPointCount;
  m_testQuadData.pointCount = testPointCount;
  cu_verify( cudaMemcpyToSymbol(constTrialQuadWeights, &trialQuadWeights[0],
      trialPointCount * sizeof(CoordinateType)) );
  cu_verify( cudaMemcpyToSymbol(constTestQuadWeights, &testQuadWeights[0],
      testPointCount * sizeof(CoordinateType)) );

  // Evaluate shapesets and copy basis data to device memory
  setupBasisData(testShapeset, trialShapeset,
      localTestQuadPoints, localTrialQuadPoints);

  if (m_isElementDataCachingEnabled == true) {

    cacheElementData(localTestQuadPoints, localTrialQuadPoints);

    if (m_isKernelDataCachingEnabled == true) {

      // Allocate memory for kernel values
      size_t size =
          maxActiveElemPairCount * trialPointCount * testPointCount;
      if (ScalarTraits<KernelType>::isComplex) size *= 2;
      m_d_kernelValues.resize(size);
    }

  } else {

    setupRawData(localTestQuadPoints, localTrialQuadPoints);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
CudaIntegrator<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    ~CudaIntegrator() {

  cu_verify( cudaSetDevice(m_deviceId) );

  thrust::device_free(m_testBasisData.values);
  thrust::device_free(m_trialBasisData.values);

  thrust::device_free(m_testBasisData.derivatives);
  thrust::device_free(m_trialBasisData.derivatives);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
void CudaIntegrator<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    integrate(
        const size_t elemPairIndexBegin, const size_t elemPairIndexEnd,
        CudaResultType *result) {

  const size_t elemPairCount = elemPairIndexEnd - elemPairIndexBegin;

  if (elemPairCount == 0)
    return;

  cu_verify( cudaSetDevice(m_deviceId) );

  // Get device pointer from host memory, no allocation or memcpy
  CudaResultType *d_result;
  cudaHostGetDevicePointer((void **)&d_result, (void *)result, 0);

  // TODO: cudaHostGetDevicePointer return non zero error code?
//  cu_verify( cudaHostGetDevicePointer((void **)&d_result, (void *)result, 0) );

  // Set kernel launch parameters
  const dim3 blockSize(m_cudaOptions.blockSize(),1,1);

  const unsigned int blockCount =
      static_cast<unsigned int>((elemPairCount-1) / blockSize.x + 1);
  const dim3 gridSize(blockCount,1,1);

  if (m_isElementDataCachingEnabled == true) {

    if (m_isKernelDataCachingEnabled == true) {

      launchCudaEvaluateIntegralFunctorKernelDataCached(
          elemPairIndexBegin, elemPairCount,
          gridSize, blockSize,
          d_result);

    } else {

      launchCudaEvaluateIntegralFunctorElementDataCached(
          elemPairIndexBegin, elemPairCount,
          gridSize, blockSize,
          d_result);
    }

  } else {

    launchCudaEvaluateIntegralFunctorNoDataCached(
        elemPairIndexBegin, elemPairCount,
        gridSize, blockSize,
        d_result);
  }
  cudaDeviceSynchronize();
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
void CudaIntegrator<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    setupBasisData(const Shapeset<BasisFunctionType> &testShapeset,
                   const Shapeset<BasisFunctionType> &trialShapeset,
                   const Matrix<CoordinateType> &localTestQuadPoints,
                   const Matrix<CoordinateType> &localTrialQuadPoints) {

  const unsigned int localCoordCount = 2;

  const unsigned int testDofCount = testShapeset.size();
  const unsigned int trialDofCount = trialShapeset.size();

  const unsigned int testPointCount = localTestQuadPoints.cols();
  const unsigned int trialPointCount = localTrialQuadPoints.cols();

  // Evaluate shapesets and copy basis data to device memory
  m_trialBasisData.dofCount = trialDofCount;
  m_testBasisData.dofCount = testDofCount;
  const unsigned int trialCompCount = 1;
  const unsigned int testCompCount = 1;
  Fiber::BasisData<BasisFunctionType> trialBasisData, testBasisData;
  size_t trialBasisDeps = 0, testBasisDeps = 0;
  trialBasisDeps |= VALUES | DERIVATIVES;
  testBasisDeps |= VALUES | DERIVATIVES;
  trialShapeset.evaluate(trialBasisDeps, localTrialQuadPoints, ALL_DOFS,
                         trialBasisData);
  testShapeset.evaluate(testBasisDeps, localTestQuadPoints, ALL_DOFS,
                        testBasisData);

  m_trialBasisData.values = thrust::device_malloc<CudaBasisFunctionType>(
      trialCompCount * trialDofCount * trialPointCount);
  m_testBasisData.values = thrust::device_malloc<CudaBasisFunctionType>(
      testCompCount * testDofCount * testPointCount);
  m_trialBasisData.derivatives = thrust::device_malloc<CudaBasisFunctionType>(
      trialCompCount * trialDofCount * trialPointCount * localCoordCount);
  m_testBasisData.derivatives = thrust::device_malloc<CudaBasisFunctionType>(
      testCompCount * testDofCount * testPointCount * localCoordCount);

  std::vector<CudaBasisFunctionType> trialBasisValues(
      trialCompCount * trialDofCount * trialPointCount);
  std::vector<CudaBasisFunctionType> testBasisValues(
      testCompCount * testDofCount * testPointCount);
  std::vector<CudaBasisFunctionType> trialBasisDerivatives(
      trialCompCount * trialDofCount * trialPointCount * localCoordCount);
  std::vector<CudaBasisFunctionType> testBasisDerivatives(
      testCompCount * testDofCount * testPointCount * localCoordCount);

  for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint)
      for (int trialComp = 0; trialComp < trialCompCount; ++trialComp)
        SetupBasisDataHelper<
        CudaBasisFunctionType, BasisFunctionType>::reorder(
            trialBasisValues, trialBasisDerivatives, trialBasisData,
            trialDofCount, trialPointCount, trialCompCount,
            trialDof, trialPoint, trialComp);
  for (int testDof = 0; testDof < testDofCount; ++testDof)
    for (int testPoint = 0; testPoint < testPointCount; ++testPoint)
      for (int testComp = 0; testComp < testCompCount; ++testComp)
        SetupBasisDataHelper<
        CudaBasisFunctionType, BasisFunctionType>::reorder(
            testBasisValues, testBasisDerivatives, testBasisData,
            testDofCount, testPointCount, testCompCount,
            testDof, testPoint, testComp);

  thrust::copy(trialBasisValues.begin(), trialBasisValues.end(),
      m_trialBasisData.values);
  thrust::copy(testBasisValues.begin(), testBasisValues.end(),
      m_testBasisData.values);
  thrust::copy(trialBasisDerivatives.begin(), trialBasisDerivatives.end(),
      m_trialBasisData.derivatives);
  thrust::copy(testBasisDerivatives.begin(), testBasisDerivatives.end(),
      m_testBasisData.derivatives);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
void CudaIntegrator<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    cacheElementData(const Matrix<CoordinateType> &localTestQuadPoints,
                     const Matrix<CoordinateType> &localTrialQuadPoints) {

  // Setup geometry data for selected elements on the device
  m_testGrid->setupGeometry();
  m_trialGrid->setupGeometry();

  // Precalculate global points on the device
  m_testGrid->local2global(
      localTestQuadPoints.template cast<CudaCoordinateType>(), m_d_testGeomData);
  m_testElemData.geomData = m_d_testGeomData.data();
  if (m_testGrid.get() == m_trialGrid.get() &&
      localTestQuadPoints == localTrialQuadPoints) {
    m_trialElemData.geomData = m_testElemData.geomData;
  } else {
    m_trialGrid->local2global(
        localTrialQuadPoints.template cast<CudaCoordinateType>(), m_d_trialGeomData);
    m_trialElemData.geomData = m_d_trialGeomData.data();
  }

  // Precalculate normals and integration elements on the device
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

  // Precalculate transposed Jacobian inverses on the device
  if ((m_kernelName == "ModifiedHelmholtz3dHypersingular" ||
      m_kernelName == "ModifiedHelmholtz3dHypersingularInterpolated")
     && fabs(m_waveNumberReal) < 1.0e-10 ) {

    m_testGrid->calculateJacobianInversesTransposed(
        m_d_testJacobianInversesTransposed);
    m_testElemData.jacobianInversesTransposed =
        m_d_testJacobianInversesTransposed.data();
    if (m_testGrid.get() == m_trialGrid.get()) {
      m_trialElemData.jacobianInversesTransposed =
          m_testElemData.jacobianInversesTransposed;
    } else {
      m_trialGrid->calculateJacobianInversesTransposed(
          m_d_trialJacobianInversesTransposed);
      m_trialElemData.jacobianInversesTransposed =
          m_d_trialJacobianInversesTransposed.data();
    }
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
void CudaIntegrator<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    setupRawData(const Matrix<CoordinateType> &localTestQuadPoints,
                 const Matrix<CoordinateType> &localTrialQuadPoints) {

  const unsigned int trialPointCount = localTrialQuadPoints.cols();
  const unsigned int testPointCount = localTestQuadPoints.cols();

  // Get raw geometry data
  m_testGrid->getRawGeometryData(
      m_testRawGeometryData.vtxCount, m_testRawGeometryData.elemCount,
      m_testRawGeometryData.vertices, m_testRawGeometryData.elementCorners);
  m_trialGrid->getRawGeometryData(
      m_trialRawGeometryData.vtxCount, m_trialRawGeometryData.elemCount,
      m_trialRawGeometryData.vertices, m_trialRawGeometryData.elementCorners);

  // Evaluate geometrical shape function values on the host
  thrust::host_vector<CoordinateType> h_testFun0(testPointCount),
      h_testFun1(testPointCount), h_testFun2(testPointCount);
  thrust::host_vector<CoordinateType> h_trialFun0(trialPointCount),
      h_trialFun1(trialPointCount), h_trialFun2(trialPointCount);

  for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
    const CoordinateType r = localTestQuadPoints(0, testPoint);
    const CoordinateType s = localTestQuadPoints(1, testPoint);
    h_testFun0[testPoint] = 1.0 - r - s;
    h_testFun1[testPoint] = r;
    h_testFun2[testPoint] = s;
  }
  for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
    const CoordinateType r = localTrialQuadPoints(0, trialPoint);
    const CoordinateType s = localTrialQuadPoints(1, trialPoint);
    h_trialFun0[trialPoint] = 1.0 - r - s;
    h_trialFun1[trialPoint] = r;
    h_trialFun2[trialPoint] = s;
  }

  // Copy geometrical shape function data to constant device memory
  cu_verify( cudaMemcpyToSymbol(constTestGeomShapeFun0, h_testFun0.data(),
      testPointCount * sizeof(CoordinateType)) );
  cu_verify( cudaMemcpyToSymbol(constTestGeomShapeFun1, h_testFun1.data(),
      testPointCount * sizeof(CoordinateType)) );
  cu_verify( cudaMemcpyToSymbol(constTestGeomShapeFun2, h_testFun2.data(),
      testPointCount * sizeof(CoordinateType)) );
  cu_verify( cudaMemcpyToSymbol(constTrialGeomShapeFun0, h_trialFun0.data(),
      trialPointCount * sizeof(CoordinateType)) );
  cu_verify( cudaMemcpyToSymbol(constTrialGeomShapeFun1, h_trialFun1.data(),
      trialPointCount * sizeof(CoordinateType)) );
  cu_verify( cudaMemcpyToSymbol(constTrialGeomShapeFun2, h_trialFun2.data(),
      trialPointCount * sizeof(CoordinateType)) );
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
void CudaIntegrator<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    launchCudaEvaluateIntegralFunctorKernelDataCached(
        const size_t elemPairIndexBegin, const size_t elemPairCount,
        const dim3 gridSize, const dim3 blockSize,
        CudaResultType *d_result) {

  // Launch kernel
  if (m_kernelName == "Laplace3dSingleLayerPotential") {

    cu_verify_void((
        CudaEvaluateLaplace3dSingleLayerPotentialKernelFunctorCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>
        <<<gridSize, blockSize>>>(
        elemPairIndexBegin, elemPairCount, m_d_testIndices.size(),
        thrust::raw_pointer_cast(m_d_testIndices.data()),
        thrust::raw_pointer_cast(m_d_trialIndices.data()),
        m_testQuadData.pointCount, m_trialQuadData.pointCount,
        m_testElemData.activeElemCount,
        thrust::raw_pointer_cast(m_testElemData.geomData),
        m_trialElemData.activeElemCount,
        thrust::raw_pointer_cast(m_trialElemData.geomData),
        thrust::raw_pointer_cast(m_d_kernelValues.data()))
    ));

  } else if (
      (m_kernelName == "ModifiedHelmholtz3dSingleLayerPotential" ||
       m_kernelName == "ModifiedHelmholtz3dSingleLayerPotentialInterpolated")
      && fabs(m_waveNumberReal) < 1.0e-10 ) {

    cu_verify_void((
        CudaEvaluateHelmholtz3dSingleLayerPotentialKernelFunctorCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>
        <<<gridSize, blockSize>>>(
        elemPairIndexBegin, elemPairCount, m_d_testIndices.size(),
        thrust::raw_pointer_cast(m_d_testIndices.data()),
        thrust::raw_pointer_cast(m_d_trialIndices.data()),
        m_testQuadData.pointCount, m_trialQuadData.pointCount,
        m_testElemData.activeElemCount,
        thrust::raw_pointer_cast(m_testElemData.geomData),
        m_trialElemData.activeElemCount,
        thrust::raw_pointer_cast(m_trialElemData.geomData),
        m_waveNumberImag, thrust::raw_pointer_cast(m_d_kernelValues.data()))
    ));

  } else {
    throw std::runtime_error(
        "CudaIntegrator::integrate(): "
        "kernel not implemented on the device");
  }
  cudaDeviceSynchronize();

  unsigned int newGridSze =
      static_cast<unsigned int>(
          (elemPairCount*m_trialBasisData.dofCount*m_testBasisData.dofCount-1)
          / blockSize.x + 1);
  dim3 newGridSize(newGridSze,1,1);

  // Launch kernel
  if (m_kernelName == "Laplace3dSingleLayerPotential") {

    cu_verify_void((
        CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorKernelCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>
        <<<newGridSize, blockSize>>>(
        elemPairIndexBegin, elemPairCount, m_d_testIndices.size(),
        thrust::raw_pointer_cast(m_d_testIndices.data()),
        thrust::raw_pointer_cast(m_d_trialIndices.data()),
        m_testQuadData.pointCount, m_trialQuadData.pointCount,
        m_testBasisData.dofCount,
        thrust::raw_pointer_cast(m_testBasisData.values),
        m_trialBasisData.dofCount,
        thrust::raw_pointer_cast(m_trialBasisData.values),
        m_testElemData.activeElemCount,
        thrust::raw_pointer_cast(m_testElemData.integrationElements),
        m_trialElemData.activeElemCount,
        thrust::raw_pointer_cast(m_trialElemData.integrationElements),
        thrust::raw_pointer_cast(m_d_kernelValues.data()),
        d_result)
    ));

  } else if (
      (m_kernelName == "ModifiedHelmholtz3dSingleLayerPotential" ||
       m_kernelName == "ModifiedHelmholtz3dSingleLayerPotentialInterpolated")
      && fabs(m_waveNumberReal) < 1.0e-10 ) {

    cu_verify_void((
        CudaEvaluateHelmholtz3dSingleLayerPotentialIntegralFunctorKernelCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>
        <<<newGridSize, blockSize>>>(
        elemPairIndexBegin, elemPairCount, m_d_testIndices.size(),
        thrust::raw_pointer_cast(m_d_testIndices.data()),
        thrust::raw_pointer_cast(m_d_trialIndices.data()),
        m_testQuadData.pointCount, m_trialQuadData.pointCount,
        m_testBasisData.dofCount,
        thrust::raw_pointer_cast(m_testBasisData.values),
        m_trialBasisData.dofCount,
        thrust::raw_pointer_cast(m_trialBasisData.values),
        m_testElemData.activeElemCount,
        thrust::raw_pointer_cast(m_testElemData.integrationElements),
        m_trialElemData.activeElemCount,
        thrust::raw_pointer_cast(m_trialElemData.integrationElements),
        thrust::raw_pointer_cast(m_d_kernelValues.data()),
        d_result)
    ));

  } else {
    throw std::runtime_error(
        "CudaIntegrator::integrate(): "
        "kernel not implemented on the device");
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
void CudaIntegrator<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    launchCudaEvaluateIntegralFunctorElementDataCached(
        const size_t elemPairIndexBegin, const size_t elemPairCount,
        const dim3 gridSize, const dim3 blockSize,
        CudaResultType *d_result) const {

  // Launch kernel
  if (m_kernelName == "Laplace3dSingleLayerPotential") {

    cu_verify_void((
        CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>
        <<<gridSize, blockSize>>>(
        elemPairIndexBegin, elemPairCount, m_d_testIndices.size(),
        thrust::raw_pointer_cast(m_d_testIndices.data()),
        thrust::raw_pointer_cast(m_d_trialIndices.data()),
        m_testQuadData.pointCount, m_trialQuadData.pointCount,
        m_testBasisData.dofCount,
        thrust::raw_pointer_cast(m_testBasisData.values),
        m_trialBasisData.dofCount,
        thrust::raw_pointer_cast(m_trialBasisData.values),
        m_testElemData.activeElemCount,
        thrust::raw_pointer_cast(m_testElemData.geomData),
        thrust::raw_pointer_cast(m_testElemData.integrationElements),
        m_trialElemData.activeElemCount,
        thrust::raw_pointer_cast(m_trialElemData.geomData),
        thrust::raw_pointer_cast(m_trialElemData.integrationElements),
        d_result)
    ));

  } else if (m_kernelName == "Laplace3dDoubleLayerPotential") {

    cu_verify_void((
        CudaEvaluateLaplace3dDoubleLayerPotentialIntegralFunctorCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>
        <<<gridSize, blockSize>>>(
        elemPairIndexBegin, elemPairCount, m_d_testIndices.size(),
        thrust::raw_pointer_cast(m_d_testIndices.data()),
        thrust::raw_pointer_cast(m_d_trialIndices.data()),
        m_testQuadData.pointCount, m_trialQuadData.pointCount,
        m_testBasisData.dofCount,
        thrust::raw_pointer_cast(m_testBasisData.values),
        m_trialBasisData.dofCount,
        thrust::raw_pointer_cast(m_trialBasisData.values),
        m_testElemData.activeElemCount,
        thrust::raw_pointer_cast(m_testElemData.geomData),
        thrust::raw_pointer_cast(m_testElemData.integrationElements),
        m_trialElemData.activeElemCount,
        thrust::raw_pointer_cast(m_trialElemData.geomData),
        thrust::raw_pointer_cast(m_trialElemData.normals),
        thrust::raw_pointer_cast(m_trialElemData.integrationElements),
        d_result)
    ));

  } else if (m_kernelName == "Laplace3dAdjointDoubleLayerPotential") {

    cu_verify_void((
        CudaEvaluateLaplace3dAdjointDoubleLayerPotentialIntegralFunctorCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>
        <<<gridSize, blockSize>>>(
        elemPairIndexBegin, elemPairCount, m_d_testIndices.size(),
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
        thrust::raw_pointer_cast(m_trialElemData.integrationElements),
        d_result)
    ));

  } else if (
      (m_kernelName == "ModifiedHelmholtz3dSingleLayerPotential" ||
       m_kernelName == "ModifiedHelmholtz3dSingleLayerPotentialInterpolated")
      && fabs(m_waveNumberReal) < 1.0e-10 ) {

    cu_verify_void((
        CudaEvaluateHelmholtz3dSingleLayerPotentialIntegralFunctorCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>
        <<<gridSize, blockSize>>>(
        elemPairIndexBegin, elemPairCount, m_d_testIndices.size(),
        thrust::raw_pointer_cast(m_d_testIndices.data()),
        thrust::raw_pointer_cast(m_d_trialIndices.data()),
        m_testQuadData.pointCount, m_trialQuadData.pointCount,
        m_testBasisData.dofCount,
        thrust::raw_pointer_cast(m_testBasisData.values),
        m_trialBasisData.dofCount,
        thrust::raw_pointer_cast(m_trialBasisData.values),
        m_testElemData.activeElemCount,
        thrust::raw_pointer_cast(m_testElemData.geomData),
        thrust::raw_pointer_cast(m_testElemData.integrationElements),
        m_trialElemData.activeElemCount,
        thrust::raw_pointer_cast(m_trialElemData.geomData),
        thrust::raw_pointer_cast(m_trialElemData.integrationElements),
        m_waveNumberImag, d_result)
    ));

  } else if (
      (m_kernelName == "ModifiedHelmholtz3dDoubleLayerPotential" ||
       m_kernelName == "ModifiedHelmholtz3dDoubleLayerPotentialInterpolated")
      && fabs(m_waveNumberReal) < 1.0e-10 ) {

    cu_verify_void((
        CudaEvaluateHelmholtz3dDoubleLayerPotentialIntegralFunctorCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>
        <<<gridSize, blockSize>>>(
        elemPairIndexBegin, elemPairCount, m_d_testIndices.size(),
        thrust::raw_pointer_cast(m_d_testIndices.data()),
        thrust::raw_pointer_cast(m_d_trialIndices.data()),
        m_testQuadData.pointCount, m_trialQuadData.pointCount,
        m_testBasisData.dofCount,
        thrust::raw_pointer_cast(m_testBasisData.values),
        m_trialBasisData.dofCount,
        thrust::raw_pointer_cast(m_trialBasisData.values),
        m_testElemData.activeElemCount,
        thrust::raw_pointer_cast(m_testElemData.geomData),
        thrust::raw_pointer_cast(m_testElemData.integrationElements),
        m_trialElemData.activeElemCount,
        thrust::raw_pointer_cast(m_trialElemData.geomData),
        thrust::raw_pointer_cast(m_trialElemData.normals),
        thrust::raw_pointer_cast(m_trialElemData.integrationElements),
        m_waveNumberImag, d_result)
    ));

  } else if (
      (m_kernelName == "ModifiedHelmholtz3dAdjointDoubleLayerPotential" ||
       m_kernelName == "ModifiedHelmholtz3dAdjointDoubleLayerPotentialInterpolated")
      && fabs(m_waveNumberReal) < 1.0e-10 ) {

    cu_verify_void((
        CudaEvaluateHelmholtz3dAdjointDoubleLayerPotentialIntegralFunctorCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>
        <<<gridSize, blockSize>>>(
        elemPairIndexBegin, elemPairCount, m_d_testIndices.size(),
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
        thrust::raw_pointer_cast(m_trialElemData.integrationElements),
        m_waveNumberImag, d_result)
    ));

  } else if (
      (m_kernelName == "ModifiedHelmholtz3dHypersingular" ||
       m_kernelName == "ModifiedHelmholtz3dHypersingularInterpolated")
      && fabs(m_waveNumberReal) < 1.0e-10 ) {

    cu_verify_void((
        CudaEvaluateHelmholtz3dHypersingularIntegralFunctorCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>
        <<<gridSize, blockSize>>>(
        elemPairIndexBegin, elemPairCount, m_d_testIndices.size(),
        thrust::raw_pointer_cast(m_d_testIndices.data()),
        thrust::raw_pointer_cast(m_d_trialIndices.data()),
        m_testQuadData.pointCount, m_trialQuadData.pointCount,
        m_testBasisData.dofCount,
        thrust::raw_pointer_cast(m_testBasisData.values),
        thrust::raw_pointer_cast(m_testBasisData.derivatives),
        m_trialBasisData.dofCount,
        thrust::raw_pointer_cast(m_trialBasisData.values),
        thrust::raw_pointer_cast(m_trialBasisData.derivatives),
        m_testElemData.activeElemCount,
        thrust::raw_pointer_cast(m_testElemData.geomData),
        thrust::raw_pointer_cast(m_testElemData.normals),
        thrust::raw_pointer_cast(m_testElemData.integrationElements),
        thrust::raw_pointer_cast(m_testElemData.jacobianInversesTransposed),
        m_trialElemData.activeElemCount,
        thrust::raw_pointer_cast(m_trialElemData.geomData),
        thrust::raw_pointer_cast(m_trialElemData.normals),
        thrust::raw_pointer_cast(m_trialElemData.integrationElements),
        thrust::raw_pointer_cast(m_trialElemData.jacobianInversesTransposed),
        m_waveNumberImag, d_result)
    ));

  } else {
    throw std::runtime_error(
        "CudaIntegrator::integrate(): "
        "kernel not implemented on the device");
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
void CudaIntegrator<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    launchCudaEvaluateIntegralFunctorNoDataCached(
        const size_t elemPairIndexBegin, const size_t elemPairCount,
        const dim3 gridSize, const dim3 blockSize,
        CudaResultType *d_result) const {

  // Launch kernel
  if (m_kernelName == "Laplace3dSingleLayerPotential") {

    cu_verify_void((
        CudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorNonCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>
        <<<gridSize, blockSize>>>(
        elemPairIndexBegin, elemPairCount, m_d_testIndices.size(),
        thrust::raw_pointer_cast(m_d_testIndices.data()),
        thrust::raw_pointer_cast(m_d_trialIndices.data()),
        m_testQuadData.pointCount, m_trialQuadData.pointCount,
        m_testBasisData.dofCount,
        thrust::raw_pointer_cast(m_testBasisData.values),
        m_trialBasisData.dofCount,
        thrust::raw_pointer_cast(m_trialBasisData.values),
        m_testRawGeometryData.elemCount, m_testRawGeometryData.vtxCount,
        thrust::raw_pointer_cast(m_testRawGeometryData.vertices),
        thrust::raw_pointer_cast(m_testRawGeometryData.elementCorners),
        m_trialRawGeometryData.elemCount, m_trialRawGeometryData.vtxCount,
        thrust::raw_pointer_cast(m_trialRawGeometryData.vertices),
        thrust::raw_pointer_cast(m_trialRawGeometryData.elementCorners),
        d_result)
    ));

  } else if (
      (m_kernelName == "ModifiedHelmholtz3dSingleLayerPotential" ||
       m_kernelName == "ModifiedHelmholtz3dSingleLayerPotentialInterpolated")
      && fabs(m_waveNumberReal) < 1.0e-10 ) {

      cu_verify_void((
          CudaEvaluateHelmholtz3dSingleLayerPotentialIntegralFunctorNonCached<
          CudaBasisFunctionType, CudaKernelType, CudaResultType>
          <<<gridSize, blockSize>>>(
          elemPairIndexBegin, elemPairCount, m_d_testIndices.size(),
          thrust::raw_pointer_cast(m_d_testIndices.data()),
          thrust::raw_pointer_cast(m_d_trialIndices.data()),
          m_testQuadData.pointCount, m_trialQuadData.pointCount,
          m_testBasisData.dofCount,
          thrust::raw_pointer_cast(m_testBasisData.values),
          m_trialBasisData.dofCount,
          thrust::raw_pointer_cast(m_trialBasisData.values),
          m_testRawGeometryData.elemCount, m_testRawGeometryData.vtxCount,
          thrust::raw_pointer_cast(m_testRawGeometryData.vertices),
          thrust::raw_pointer_cast(m_testRawGeometryData.elementCorners),
          m_trialRawGeometryData.elemCount, m_trialRawGeometryData.vtxCount,
          thrust::raw_pointer_cast(m_trialRawGeometryData.vertices),
          thrust::raw_pointer_cast(m_trialRawGeometryData.elementCorners),
          m_waveNumberImag, d_result)
      ));

  } else {
    throw std::runtime_error(
        "CudaIntegrator::integrate(): "
        "kernel not implemented on the device");
  }
}

// Explicit instantiations
template class CudaIntegrator<float, float, float, float, float, float>;
template class CudaIntegrator<float, float, std::complex<float>, float, float, float>;
template class CudaIntegrator<float, std::complex<float>, std::complex<float>, float, float, float>;
template class CudaIntegrator<std::complex<float>, float, std::complex<float>, float, float, float>;
template class CudaIntegrator<std::complex<float>, std::complex<float>, std::complex<float>, float, float, float>;

template class CudaIntegrator<double, double, double, double, double, double>;
template class CudaIntegrator<double, double, std::complex<double>, double, double, double>;
template class CudaIntegrator<double, std::complex<double>, std::complex<double>, double, double, double>;
template class CudaIntegrator<std::complex<double>, double, std::complex<double>, double, double, double>;
template class CudaIntegrator<std::complex<double>, std::complex<double>, std::complex<double>, double, double, double>;

template class CudaIntegrator<float, float, float, double, double, double>;
template class CudaIntegrator<float, float, std::complex<float>, double, double, double>;
template class CudaIntegrator<float, std::complex<float>, std::complex<float>, double, double, double>;
template class CudaIntegrator<std::complex<float>, float, std::complex<float>, double, double, double>;
template class CudaIntegrator<std::complex<float>, std::complex<float>, std::complex<float>, double, double, double>;

template class CudaIntegrator<double, double, double, float, float, float>;
template class CudaIntegrator<double, double, std::complex<double>, float, float, float>;
template class CudaIntegrator<double, std::complex<double>, std::complex<double>, float, float, float>;
template class CudaIntegrator<std::complex<double>, double, std::complex<double>, float, float, float>;
template class CudaIntegrator<std::complex<double>, std::complex<double>, std::complex<double>, float, float, float>;

} // namespace Fiber
