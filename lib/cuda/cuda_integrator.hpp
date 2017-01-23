// Copyright (C) 2011-2012 by the Bem++ Authors
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

#ifndef fiber_cuda_integrator_hpp
#define fiber_cuda_integrator_hpp

#include "cuda_options.hpp"
#include "cuda.hpp"

#include "../common/common.hpp"
#include "../common/eigen_support.hpp"
#include "../common/scalar_traits.hpp"

#include "../fiber/types.hpp"
#include "../fiber/shared_ptr.hpp"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

namespace Bempp {
/** \cond FORWARD_DECL */
template <typename CoordinateType> class CudaGrid;
/** \endcond */

} // namespace Bempp

namespace Fiber {
/** \cond FORWARD_DECL */
template <typename BasisFunctionType> class Shapeset;
template <typename BasisFunctionType> class BasisData;
template <typename KernelType> class CollectionOfKernels;
/** \endcond */

/** \brief Regular integration over pairs of elements on the device. */
template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
class CudaIntegrator {
public:

  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
  typedef typename ScalarTraits<CudaResultType>::RealType CudaCoordinateType;

  /** \brief Constructor */
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
      const int deviceId, const Bempp::CudaOptions &cudaOptions);

  CudaIntegrator(
      const Matrix<CoordinateType> &localTestQuadPoints,
      const Matrix<CoordinateType> &localTrialQuadPoints,
      const std::vector<CoordinateType> &testQuadWeights,
      const std::vector<CoordinateType> &trialQuadWeights,
      const Shapeset<BasisFunctionType> &testShapeset,
      const Shapeset<BasisFunctionType> &trialShapeset,
      const shared_ptr<Bempp::CudaGrid<CudaCoordinateType>> &testGrid,
      const shared_ptr<Bempp::CudaGrid<CudaCoordinateType>> &trialGrid,
      const size_t maxActiveElemPairCount,
      const shared_ptr<const CollectionOfKernels<KernelType>> &kernel,
      const int deviceId, const Bempp::CudaOptions &cudaOptions);

  /** \brief Destructor. */
  virtual ~CudaIntegrator();

  void integrate(const size_t elemPairIndexBegin,
                 const size_t elemPairIndexEnd,
                 CudaResultType *result);

  void integrate(const size_t elemPairCount,
                 CudaResultType *result);

  void pushElemPairIndicesToDevice(
      const std::vector<int> &testElemPairIndices,
      const std::vector<int> &trialElemPairIndices,
      const size_t elemPairCount);

private:
  /** \cond PRIVATE */
  void setupBasisData(const Shapeset<BasisFunctionType> &testShapeset,
                      const Shapeset<BasisFunctionType> &trialShapeset,
                      const Matrix<CoordinateType> &localTestQuadPoints,
                      const Matrix<CoordinateType> &localTrialQuadPoints);

  void cacheElementData(const Matrix<CoordinateType> &localTestQuadPoints,
                        const Matrix<CoordinateType> &localTrialQuadPoints);

  void setupRawData(const Matrix<CoordinateType> &localTestQuadPoints,
                    const Matrix<CoordinateType> &localTrialQuadPoints);

  void launchCudaEvaluateIntegralFunctorKernelDataCached(
      const size_t elemPairIndexBegin, const size_t elemPairCount,
      const dim3 gridSize, const dim3 blockSize,
      CudaResultType *d_result);

  void launchCudaEvaluateIntegralFunctorKernelDataCached(
      const size_t elemPairCount,
      const dim3 gridSize, const dim3 blockSize,
      CudaResultType *d_result);

  void launchCudaEvaluateIntegralFunctorElementDataCached(
      const size_t elemPairIndexBegin, const size_t elemPairCount,
      const dim3 gridSize, const dim3 blockSize,
      CudaResultType *d_result) const;

  void launchCudaEvaluateIntegralFunctorElementDataCached(
      const size_t elemPairCount,
      const dim3 gridSize, const dim3 blockSize,
      CudaResultType *d_result) const;

  void launchCudaEvaluateIntegralFunctorNoDataCached(
      const size_t elemPairIndexBegin, const size_t elemPairCount,
      const dim3 gridSize, const dim3 blockSize,
      CudaResultType *d_result) const;

  void launchCudaEvaluateIntegralFunctorNoDataCached(
      const size_t elemPairCount,
      const dim3 gridSize, const dim3 blockSize,
      CudaResultType *d_result) const;

  const int m_deviceId;
  const std::string m_kernelName;
  const Bempp::CudaOptions m_cudaOptions;
  const bool m_isElementDataCachingEnabled, m_isKernelDataCachingEnabled;
  const shared_ptr<Bempp::CudaGrid<CudaCoordinateType>> m_testGrid, m_trialGrid;

  CudaKernelType m_waveNumberReal, m_waveNumberImag;

  thrust::device_vector<int> m_d_testIndices, m_d_trialIndices;

  QuadData<CudaCoordinateType> m_testQuadData, m_trialQuadData;
  BasisFunData<CudaBasisFunctionType> m_testBasisData, m_trialBasisData;
  ElemData<CudaBasisFunctionType> m_testElemData, m_trialElemData;
  thrust::device_vector<CudaCoordinateType> m_d_testGeomData, m_d_trialGeomData;
  thrust::device_vector<CudaCoordinateType> m_d_testNormals, m_d_trialNormals;
  thrust::device_vector<CudaCoordinateType>
  m_d_testIntegrationElements, m_d_trialIntegrationElements;
  thrust::device_vector<CudaCoordinateType>
  m_d_testSurfaceCurls, m_d_trialSurfaceCurls;
  RawGeometryData<CudaCoordinateType>
  m_testRawGeometryData, m_trialRawGeometryData;

  thrust::device_vector<CudaKernelType> m_d_kernelValues;

  /** \endcond */
};

} // namespace Fiber

#endif
