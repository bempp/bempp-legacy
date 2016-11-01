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

namespace Bempp {

/** \cond FORWARD_DECL */
class CudaGrid;
/** \endcond */

} // namespace Bempp

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename BasisFunctionType> class Shapeset;
template <typename BasisFunctionType> class BasisData;
/** \endcond */

/** \brief Regular integration over pairs of elements on the device. */
template <typename BasisFunctionType, typename ResultType>
class CudaIntegrator {
public:

  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

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
      bool cacheElemData = false);

  /** \brief Destructor. */
  virtual ~CudaIntegrator();

  void integrate(
      const std::vector<int> &elementPairTestIndices,
      const std::vector<int> &elementPairTrialIndices,
      std::vector<Matrix<ResultType>*> &result);

private:
  /** \cond PRIVATE */

  QuadData<CoordinateType> m_testQuadData;
  QuadData<CoordinateType> m_trialQuadData;

  BasisFunData<BasisFunctionType> m_testBasisData;
  BasisFunData<BasisFunctionType> m_trialBasisData;

  Matrix<CoordinateType> m_localTestQuadPoints;
  Matrix<CoordinateType> m_localTrialQuadPoints;

  shared_ptr<Bempp::CudaGrid> m_testGrid;
  shared_ptr<Bempp::CudaGrid> m_trialGrid;

  bool m_cacheElemData;

  /** \endcond */
};

} // namespace Fiber

#endif
