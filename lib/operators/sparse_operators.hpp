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

#ifndef bempp_sparse_operators_hpp
#define bempp_sparse_operators_hpp

#include "../assembly/symmetry.hpp"
#include "../common/types.hpp"
#include "../common/shared_ptr.hpp"
#include "../assembly/elementary_local_operator.hpp"
#include "../assembly/general_elementary_local_operator.hpp"
#include "../operators/abstract_identity_operator.hpp"
#include "../operators/abstract_maxwell_identity_operator.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../fiber/default_test_trial_integral_imp.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/simple_test_trial_integrand_functor.hpp"
#include "../fiber/surface_grad_3d_functor.hpp"
#include "../fiber/shared_ptr.hpp"
#include "../assembly/context.hpp"

#include <boost/type_traits/is_complex.hpp>

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
class ElementaryLocalOperator;
template <typename BasisFunctionType> class Space;

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const ElementaryLocalOperator<BasisFunctionType, ResultType>>
identityOperator(const ParameterList &parameterList,
                 const shared_ptr<const Space<BasisFunctionType>> &domain,
                 const shared_ptr<const Space<BasisFunctionType>> &range,
                 const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
                 const std::string &label = "", int symmetry = NO_SYMMETRY) {

  return shared_ptr<
      const ElementaryLocalOperator<BasisFunctionType, ResultType>>(
      new AbstractIdentityOperator<BasisFunctionType, ResultType>(
          domain, range, dualToRange, label, symmetry));
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const ElementaryLocalOperator<BasisFunctionType, ResultType>>
maxwellIdentityOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label = "", int symmetry = NO_SYMMETRY) {

  return shared_ptr<
      const ElementaryLocalOperator<BasisFunctionType, ResultType>>(
      new AbstractMaxwellIdentityOperator<BasisFunctionType, ResultType>(
          domain, range, dualToRange, label, symmetry));
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const ElementaryLocalOperator<BasisFunctionType, ResultType>>
laplaceBeltramiOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label = "", int symmetry = NO_SYMMETRY) {

  if (domain->codomainDimension() != 1)
    throw std::invalid_argument(
        "laplaceBeltrami3dOperator2(): "
        "domain must consist of scalar-valued functions");
  if (range->codomainDimension() != 1)
    throw std::invalid_argument(
        "laplaceBeltrami3dOperator2(): "
        "range must consist of scalar-valued functions");
  if (dualToRange->codomainDimension() != 1)
    throw std::invalid_argument(
        "laplaceBeltrami3dOperator2(): "
        "space dual to range must consist of scalar-valued functions");

  typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

  typedef Fiber::SurfaceGrad3dFunctor<CoordinateType> TransformationFunctor;
  typedef Fiber::SimpleTestTrialIntegrandFunctor<BasisFunctionType, ResultType>
      IntegrandFunctor;

  typedef GeneralElementaryLocalOperator<BasisFunctionType, ResultType> Op;
  return shared_ptr<
      const ElementaryLocalOperator<BasisFunctionType, ResultType>>(
      new GeneralElementaryLocalOperator<BasisFunctionType, ResultType>(
          domain, range, dualToRange, label, symmetry, TransformationFunctor(),
          TransformationFunctor(), IntegrandFunctor()));
}
}

#endif
