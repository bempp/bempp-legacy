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
#include "cuda_laplace_3d_single_layer_potential_kernel_functor.hpp"
#include "cuda_laplace_3d_double_layer_potential_kernel_functor.hpp"
#include "cuda_laplace_3d_adjoint_double_layer_potential_kernel_functor.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/shapeset.hpp"
#include "../fiber/basis_data.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../fiber/collection_of_kernels.hpp"

namespace Fiber {

template <typename BasisFunctionType, typename KernelType, typename ResultType,
    typename KernelFunctor>
CudaIntegrator<BasisFunctionType, KernelType, ResultType, KernelFunctor>::
    CudaIntegrator(
        const Matrix<CoordinateType> &localTestQuadPoints,
        const Matrix<CoordinateType> &localTrialQuadPoints,
        const std::vector<CoordinateType> &testQuadWeights,
        const std::vector<CoordinateType> &trialQuadWeights,
        const Shapeset<BasisFunctionType> &testShapeset,
        const Shapeset<BasisFunctionType> &trialShapeset,
        shared_ptr<Bempp::CudaGrid<CoordinateType>> testGrid,
        shared_ptr<Bempp::CudaGrid<CoordinateType>> trialGrid,
        const std::vector<int> &testIndices, const std::vector<int> &trialIndices,
        const shared_ptr<const CollectionOfKernels<KernelType>> &kernel,
        const int deviceId, const Bempp::CudaOptions &cudaOptions)
        : m_testGrid(testGrid), m_trialGrid(trialGrid),
          m_deviceId(deviceId), m_cudaOptions(cudaOptions) {

  const size_t trialIndexCount = trialIndices.size();
  const size_t testIndexCount = testIndices.size();

  const unsigned int trialDofCount = trialShapeset.size();
  const unsigned int testDofCount = testShapeset.size();

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

  cu_verify( cudaSetDevice(m_deviceId) );

  // Get device function pointer for kernel
  cu_verify( cudaMalloc((void **)&m_d_kernel, sizeof(KernelFunctor)) );
  cu_verify( cudaMemcpy(m_d_kernel, &(*(kernel->cudaFunctor())),
      sizeof(KernelFunctor), cudaMemcpyHostToDevice) );
  kernel->addGeometricalDependencies(m_testGeomDeps, m_trialGeomDeps);

  // Copy element indices to device memory
  m_d_trialIndices.resize(trialIndexCount);
  m_d_testIndices.resize(testIndexCount);
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
  m_trialBasisData.dofCount = trialDofCount;
  m_testBasisData.dofCount = testDofCount;
  const unsigned int trialCompCount = 1;
  const unsigned int testCompCount = 1;
  BasisData<BasisFunctionType> trialBasisData, testBasisData;
  size_t trialBasisDeps = 0, testBasisDeps = 0;
  trialBasisDeps |= VALUES;
  testBasisDeps |= VALUES;
  trialShapeset.evaluate(trialBasisDeps, localTrialQuadPoints, ALL_DOFS,
                         trialBasisData);
  testShapeset.evaluate(testBasisDeps, localTestQuadPoints, ALL_DOFS,
                        testBasisData);
  m_trialBasisData.values = thrust::device_malloc<CudaBasisFunctionType>(
      trialCompCount * trialDofCount * trialPointCount);
  m_testBasisData.values = thrust::device_malloc<CudaBasisFunctionType>(
      testCompCount * testDofCount * testPointCount);
//  thrust::copy(trialBasisData.values.begin(), trialBasisData.values.end(),
//      m_trialBasisData.values);
//  thrust::copy(testBasisData.values.begin(), testBasisData.values.end(),
//      m_testBasisData.values);
  std::vector<CudaBasisFunctionType> trialBasisValues(
      trialCompCount * trialDofCount * trialPointCount);
  std::vector<CudaBasisFunctionType> testBasisValues(
      testCompCount * testDofCount * testPointCount);
  for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint)
      for (int trialComp = 0; trialComp < trialCompCount; ++trialComp)
        trialBasisValues[trialDof * trialPointCount * trialCompCount
                        + trialPoint * trialCompCount
                        + trialComp] =
            trialBasisData.values(trialComp, trialDof, trialPoint);
  for (int testDof = 0; testDof < testDofCount; ++testDof)
    for (int testPoint = 0; testPoint < testPointCount; ++testPoint)
      for (int testComp = 0; testComp < testCompCount; ++testComp)
        testBasisValues[testDof * testPointCount * testCompCount
                        + testPoint * testCompCount
                        + testComp] =
            testBasisData.values(testComp, testDof, testPoint);
  thrust::copy(trialBasisValues.begin(), trialBasisValues.end(),
      m_trialBasisData.values);
  thrust::copy(testBasisValues.begin(), testBasisValues.end(),
      m_testBasisData.values);
//  cu_verify( cudaMemcpyToSymbol(constTrialBasisValues, &(*(trialBasisValues.begin())),
//      trialBasisValues.size() * sizeof(CudaBasisFunctionType)) );
//  cu_verify( cudaMemcpyToSymbol(constTestBasisValues, &(*(testBasisValues.begin())),
//      testBasisValues.size() * sizeof(CudaBasisFunctionType)) );

  if (cudaOptions.isElementDataCachingEnabled() == true) {

    // Setup geometry data for selected elements on the device
    m_testGrid->setupGeometry();
    m_trialGrid->setupGeometry();

    // Precalculate global points on the device
    m_testGrid->local2global(localTestQuadPoints, m_d_testGeomData);
    m_testElemData.geomData = m_d_testGeomData.data();
    if (m_testGrid.get() == m_trialGrid.get() &&
        localTestQuadPoints == localTrialQuadPoints) {
      m_trialElemData.geomData = m_testElemData.geomData;
    } else {
      m_trialGrid->local2global(localTrialQuadPoints, m_d_trialGeomData);
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

  } else {

    // Get raw geometry data
    m_testGrid->getRawGeometryData(
        m_testRawGeometryData.vtxCount, m_testRawGeometryData.elemCount,
        m_testRawGeometryData.vertices, m_testRawGeometryData.elementCorners);
    m_trialGrid->getRawGeometryData(
        m_trialRawGeometryData.vtxCount, m_trialRawGeometryData.elemCount,
        m_trialRawGeometryData.vertices, m_trialRawGeometryData.elementCorners);

    // Evaluate geometrical shape function values on the host
    thrust::host_vector<CoordinateType> h_testFun0(testPointCount);
    thrust::host_vector<CoordinateType> h_testFun1(testPointCount);
    thrust::host_vector<CoordinateType> h_testFun2(testPointCount);
    thrust::host_vector<CoordinateType> h_trialFun0(trialPointCount);
    thrust::host_vector<CoordinateType> h_trialFun1(trialPointCount);
    thrust::host_vector<CoordinateType> h_trialFun2(trialPointCount);

    for (int testPoint = 0; testPoint < testPointCount; ++testPoint) {
      const CoordinateType r = localTestQuadPoints(0,testPoint);
      const CoordinateType s = localTestQuadPoints(1,testPoint);
      h_testFun0[testPoint] = 1.0 - r - s;
      h_testFun1[testPoint] = r;
      h_testFun2[testPoint] = s;
    }
    for (int trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
      const CoordinateType r = localTrialQuadPoints(0,trialPoint);
      const CoordinateType s = localTrialQuadPoints(1,trialPoint);
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
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
    typename KernelFunctor>
CudaIntegrator<BasisFunctionType, KernelType, ResultType, KernelFunctor>::
    ~CudaIntegrator() {

  cu_verify( cudaSetDevice(m_deviceId) );

  cudaFree(m_d_kernel);

  thrust::device_free(m_testBasisData.values);
  thrust::device_free(m_trialBasisData.values);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
    typename KernelFunctor>
void CudaIntegrator<BasisFunctionType, KernelType, ResultType, KernelFunctor>::
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
  unsigned int blockSze = m_cudaOptions.blockSize();
  dim3 blockSize(blockSze,1,1);

  unsigned int gridSze =
      static_cast<unsigned int>((elemPairCount-1) / blockSze + 1);
  dim3 gridSize(gridSze,1,1);

  if (m_cudaOptions.isElementDataCachingEnabled() == true) {

    if (m_cudaOptions.isKernelDataCachingEnabled() == true) {

      thrust::device_vector<CudaKernelType> d_kernelValues(
          elemPairCount * m_trialQuadData.pointCount * m_testQuadData.pointCount);

//      unsigned int newerGridSze =
//          static_cast<unsigned int>(
//              (elemPairCount*m_testQuadData.pointCount*m_trialQuadData.pointCount-1)
//              / blockSze + 1);
//      dim3 newerGridSize(newerGridSze,1,1);

      // Launch kernel
      cu_verify_void((
          RawCudaEvaluateLaplace3dSingleLayerPotentialKernelFunctorCached<
          CudaBasisFunctionType, CudaKernelType, CudaResultType>
          <<<gridSize, blockSize>>>(
          elemPairIndexBegin, elemPairCount, m_d_testIndices.size(),
          thrust::raw_pointer_cast(m_d_testIndices.data()),
          thrust::raw_pointer_cast(m_d_trialIndices.data()),
          m_testQuadData.pointCount, m_trialQuadData.pointCount,
          m_testElemData.activeElemCount,
          thrust::raw_pointer_cast(m_testElemData.geomData),
          thrust::raw_pointer_cast(m_testElemData.normals),
          m_trialElemData.activeElemCount,
          thrust::raw_pointer_cast(m_trialElemData.geomData),
          thrust::raw_pointer_cast(m_trialElemData.normals),
          thrust::raw_pointer_cast(d_kernelValues.data()),
          *m_d_kernel, m_testGeomDeps, m_trialGeomDeps)
      ));
      cudaDeviceSynchronize();

      unsigned int newGridSze =
          static_cast<unsigned int>(
              (elemPairCount*m_trialBasisData.dofCount*m_testBasisData.dofCount-1)
              / blockSze + 1);
      dim3 newGridSize(newGridSze,1,1);

      // Launch kernel
      cu_verify_void((
          RawCudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorKernelCached<
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
          thrust::raw_pointer_cast(d_kernelValues.data()),
          d_result)
      ));

    } else {

      // Launch kernel
      cu_verify_void((
          RawCudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorCached<
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
          thrust::raw_pointer_cast(m_trialElemData.normals),
          thrust::raw_pointer_cast(m_trialElemData.integrationElements),
          *m_d_kernel, m_testGeomDeps, m_trialGeomDeps,
          d_result)
      ));
    }

  } else {

    // Launch kernel
    cu_verify_void((
        RawCudaEvaluateLaplace3dSingleLayerPotentialIntegralFunctorNonCached<
        CudaBasisFunctionType, CudaKernelType, CudaResultType>
        <<<gridSize, blockSize>>>(
        elemPairIndexBegin, elemPairCount, m_d_trialIndices.size(),
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
        *m_d_kernel, m_testGeomDeps, m_trialGeomDeps,
        d_result)
    ));
  }
  cudaDeviceSynchronize();
}

// Explicit instantiations
template <typename KernelType> class CudaLaplace3dSingleLayerPotentialKernelFunctor;
template class CudaIntegrator<double, double, double, CudaLaplace3dSingleLayerPotentialKernelFunctor<double>>;
template class CudaIntegrator<double, std::complex<double>, std::complex<double>, CudaLaplace3dSingleLayerPotentialKernelFunctor<thrust::complex<double>>>;
template class CudaIntegrator<std::complex<double>, std::complex<double>, std::complex<double>, CudaLaplace3dSingleLayerPotentialKernelFunctor<thrust::complex<double>>>;
template class CudaIntegrator<float, float, float, CudaLaplace3dSingleLayerPotentialKernelFunctor<float>>;
template class CudaIntegrator<float, std::complex<float>, std::complex<float>, CudaLaplace3dSingleLayerPotentialKernelFunctor<thrust::complex<float>>>;
template class CudaIntegrator<std::complex<float>, std::complex<float>, std::complex<float>, CudaLaplace3dSingleLayerPotentialKernelFunctor<thrust::complex<float>>>;

template <typename KernelType> class CudaLaplace3dDoubleLayerPotentialKernelFunctor;
template class CudaIntegrator<double, double, double, CudaLaplace3dDoubleLayerPotentialKernelFunctor<double>>;
//template class CudaIntegrator<double, double, std::complex<double>, CudaLaplace3dDoubleLayerPotentialKernelFunctor<double>>;
//template class CudaIntegrator<double, std::complex<double>, std::complex<double>, CudaLaplace3dDoubleLayerPotentialKernelFunctor<thrust::complex<double>>>;
//template class CudaIntegrator<std::complex<double>, double, std::complex<double>, CudaLaplace3dDoubleLayerPotentialKernelFunctor<double>>;
//template class CudaIntegrator<std::complex<double>, std::complex<double>, std::complex<double>, CudaLaplace3dDoubleLayerPotentialKernelFunctor<thrust::complex<double>>>;
template class CudaIntegrator<float, float, float, CudaLaplace3dDoubleLayerPotentialKernelFunctor<float>>;
//template class CudaIntegrator<float, float, std::complex<float>, CudaLaplace3dDoubleLayerPotentialKernelFunctor<float>>;
//template class CudaIntegrator<float, std::complex<float>, std::complex<float>, CudaLaplace3dDoubleLayerPotentialKernelFunctor<thrust::complex<float>>>;
//template class CudaIntegrator<std::complex<float>, float, std::complex<float>, CudaLaplace3dDoubleLayerPotentialKernelFunctor<float>>;
//template class CudaIntegrator<std::complex<float>, std::complex<float>, std::complex<float>, CudaLaplace3dDoubleLayerPotentialKernelFunctor<thrust::complex<float>>>;

template <typename KernelType> class CudaLaplace3dAdjointDoubleLayerPotentialKernelFunctor;
template class CudaIntegrator<double, double, double, CudaLaplace3dAdjointDoubleLayerPotentialKernelFunctor<double>>;
//template class CudaIntegrator<double, double, std::complex<double>, CudaLaplace3dAdjointDoubleLayerPotentialKernelFunctor<double>>;
//template class CudaIntegrator<double, std::complex<double>, std::complex<double>, CudaLaplace3dAdjointDoubleLayerPotentialKernelFunctor<thrust::complex<double>>>;
//template class CudaIntegrator<std::complex<double>, double, std::complex<double>, CudaLaplace3dAdjointDoubleLayerPotentialKernelFunctor<double>>;
//template class CudaIntegrator<std::complex<double>, std::complex<double>, std::complex<double>, CudaLaplace3dAdjointDoubleLayerPotentialKernelFunctor<thrust::complex<double>>>;
template class CudaIntegrator<float, float, float, CudaLaplace3dAdjointDoubleLayerPotentialKernelFunctor<float>>;
//template class CudaIntegrator<float, float, std::complex<float>, CudaLaplace3dAdjointDoubleLayerPotentialKernelFunctor<float>>;
//template class CudaIntegrator<float, std::complex<float>, std::complex<float>, CudaLaplace3dAdjointDoubleLayerPotentialKernelFunctor<thrust::complex<float>>>;
//template class CudaIntegrator<std::complex<float>, float, std::complex<float>, CudaLaplace3dAdjointDoubleLayerPotentialKernelFunctor<float>>;
//template class CudaIntegrator<std::complex<float>, std::complex<float>, std::complex<float>, CudaLaplace3dAdjointDoubleLayerPotentialKernelFunctor<thrust::complex<float>>>;

} // namespace Fiber
