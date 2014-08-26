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

#include "identity_operator.hpp"

#include "boundary_operator.hpp"

#include "../common/boost_make_shared_fwd.hpp"
#include "../fiber/default_test_trial_integral_imp.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/simple_test_trial_integrand_functor.hpp"

#include <boost/type_traits/is_complex.hpp>

namespace Bempp {

////////////////////////////////////////////////////////////////////////////////
// IdentityOperatorId

template <typename BasisFunctionType, typename ResultType>
IdentityOperatorId<BasisFunctionType, ResultType>::IdentityOperatorId(
    const IdentityOperator<BasisFunctionType, ResultType> &op)
    : m_domain(op.domain().get()), m_range(op.range().get()),
      m_dualToRange(op.dualToRange().get()) {}

template <typename BasisFunctionType, typename ResultType>
size_t IdentityOperatorId<BasisFunctionType, ResultType>::hash() const {
  typedef IdentityOperator<BasisFunctionType, ResultType> OperatorType;
  size_t result = tbb::tbb_hasher(typeid(OperatorType).name());
  tbb_hash_combine(result, m_domain);
  tbb_hash_combine(result, m_range);
  tbb_hash_combine(result, m_dualToRange);
  return result;
}

template <typename BasisFunctionType, typename ResultType>
bool IdentityOperatorId<BasisFunctionType, ResultType>::isEqual(
    const AbstractBoundaryOperatorId &other) const {
  // dynamic_cast won't suffice since we want to make sure both objects
  // are of exactly the same type (dynamic_cast would succeed for a subclass)
  if (typeid(other) == typeid(*this)) {
    const IdentityOperatorId &otherCompatible =
        static_cast<const IdentityOperatorId &>(other);
    return (m_domain == otherCompatible.m_domain &&
            m_range == otherCompatible.m_range &&
            m_dualToRange == otherCompatible.m_dualToRange);
  } else
    return false;
}

////////////////////////////////////////////////////////////////////////////////
// IdentityOperator

template <typename BasisFunctionType, typename ResultType>
IdentityOperator<BasisFunctionType, ResultType>::IdentityOperator(
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label, int symmetry)
    : Base(domain, range, dualToRange, label,
           (symmetry & AUTO_SYMMETRY
                ? (domain == dualToRange
                       ? (boost::is_complex<BasisFunctionType>()
                              ? HERMITIAN
                              : SYMMETRIC | HERMITIAN)
                       : NO_SYMMETRY)
                : symmetry)),
      m_id(
          boost::make_shared<IdentityOperatorId<BasisFunctionType, ResultType>>(
              *this)) {
  typedef Fiber::SimpleTestTrialIntegrandFunctor<BasisFunctionType, ResultType>
  IntegrandFunctor;
  typedef Fiber::DefaultTestTrialIntegral<IntegrandFunctor> Integral;
  m_integral.reset(new Integral((IntegrandFunctor())));
}

template <typename BasisFunctionType, typename ResultType>
IdentityOperator<BasisFunctionType, ResultType>::~IdentityOperator() {}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const AbstractBoundaryOperatorId>
IdentityOperator<BasisFunctionType, ResultType>::id() const {
  return m_id;
}

template <typename BasisFunctionType, typename ResultType>
const typename IdentityOperator<
    BasisFunctionType, ResultType>::CollectionOfShapesetTransformations &
IdentityOperator<BasisFunctionType, ResultType>::testTransformations() const {
  return this->dualToRange()->basisFunctionValue();
}

template <typename BasisFunctionType, typename ResultType>
const typename IdentityOperator<
    BasisFunctionType, ResultType>::CollectionOfShapesetTransformations &
IdentityOperator<BasisFunctionType, ResultType>::trialTransformations() const {
  return this->domain()->basisFunctionValue();
}

template <typename BasisFunctionType, typename ResultType>
const typename IdentityOperator<BasisFunctionType,
                                ResultType>::TestTrialIntegral &
IdentityOperator<BasisFunctionType, ResultType>::integral() const {
  return *m_integral;
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType> identityOperator(
    const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label, int symmetry) {
  typedef IdentityOperator<BasisFunctionType, ResultType> Id;
  return BoundaryOperator<BasisFunctionType, ResultType>(
      context,
      boost::make_shared<Id>(domain, range, dualToRange, label, symmetry));
}

#define INSTANTIATE_NONMEMBER_CONSTRUCTOR(BASIS, RESULT)                       \
  template BoundaryOperator<BASIS, RESULT> identityOperator(                   \
      const shared_ptr<const Context<BASIS, RESULT>> &,                        \
      const shared_ptr<const Space<BASIS>> &,                                  \
      const shared_ptr<const Space<BASIS>> &,                                  \
      const shared_ptr<const Space<BASIS>> &, const std::string &, int)
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(IdentityOperator);

} // namespace Bempp
