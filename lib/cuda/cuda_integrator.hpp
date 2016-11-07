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

#include "cuda_evaluate_laplace_3d_single_layer_potential_integral_functor.cuh"
#include "cuda_evaluate_laplace_3d_double_layer_potential_integral_functor.cuh"

#include "../common/common.hpp"
#include "../common/eigen_support.hpp"
#include "../common/scalar_traits.hpp"

#include "../fiber/types.hpp"
#include "../fiber/shared_ptr.hpp"

#include <thrust/device_ptr.h>
#include <thrust/host_vector.h>
#include <thrust/complex.h>

#include <type_traits>

namespace Bempp {

/** \cond FORWARD_DECL */
class CudaGrid;
/** \endcond */

} // namespace Bempp

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename BasisFunctionType> class Shapeset;
template <typename BasisFunctionType> class BasisData;
template <typename KernelType> class CollectionOfKernels;
/** \endcond */

/** \brief Regular integration over pairs of elements on the device. */
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class CudaIntegrator {
public:

  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  // TODO: Check if this is correct
  typedef typename thrust::complex<CoordinateType> CudaComplexType;
  typedef typename std::conditional<
      std::is_same<BasisFunctionType,CoordinateType>::value,
      BasisFunctionType, CudaComplexType>::type
      CudaBasisFunctionType;
  typedef typename std::conditional<
      std::is_same<KernelType,CoordinateType>::value,
      KernelType, CudaComplexType>::type
      CudaKernelType;
  typedef typename std::conditional<
      std::is_same<ResultType,CoordinateType>::value,
      ResultType, CudaComplexType>::type
      CudaResultType;

  /** \brief Constructor */
  CudaIntegrator(
      const Matrix<CoordinateType> &localTestQuadPoints,
      const Matrix<CoordinateType> &localTrialQuadPoints,
      const std::vector<CoordinateType> &testQuadWeights,
      const std::vector<CoordinateType> &trialQuadWeights,
      const Shapeset<BasisFunctionType> &testShapeset,
      const Shapeset<BasisFunctionType> &trialShapeset,
      shared_ptr<Bempp::CudaGrid> testGrid,
      shared_ptr<Bempp::CudaGrid> trialGrid,
      const CollectionOfKernels<KernelType> &kernels,
      bool cacheElemData = false);

  /** \brief Destructor. */
  virtual ~CudaIntegrator();

  void integrate(
      std::vector<int>::iterator startElementPairTestIndices,
      std::vector<int>::iterator endElementPairTestIndices,
      std::vector<int>::iterator startElementPairTrialIndices,
      std::vector<int>::iterator endElementPairTrialIndices,
      typename thrust::host_vector<ResultType>::iterator startResult,
      typename thrust::host_vector<ResultType>::iterator endResult);

private:
  /** \cond PRIVATE */

  QuadData<CoordinateType> m_testQuadData;
  QuadData<CoordinateType> m_trialQuadData;

  BasisFunData<CudaBasisFunctionType> m_testBasisData;
  BasisFunData<CudaBasisFunctionType> m_trialBasisData;

  Matrix<CoordinateType> m_localTestQuadPoints;
  Matrix<CoordinateType> m_localTrialQuadPoints;

  shared_ptr<Bempp::CudaGrid> m_testGrid;
  shared_ptr<Bempp::CudaGrid> m_trialGrid;

  const CollectionOfKernels<KernelType> &m_kernels;

  bool m_cacheElemData;

  /** \endcond */
};

} // namespace Fiber

#endif
