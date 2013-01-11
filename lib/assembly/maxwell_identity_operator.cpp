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

#include "maxwell_identity_operator.hpp"

#include "boundary_operator.hpp"

#include "../common/boost_make_shared_fwd.hpp"
#include "../fiber/default_test_trial_integral_imp.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/maxwell_test_trial_integrand_functor.hpp"

#include <boost/type_traits/is_complex.hpp>

namespace Bempp
{

////////////////////////////////////////////////////////////////////////////////
// MaxwellIdentityOperator

template <typename BasisFunctionType, typename ResultType>
MaxwellIdentityOperator<BasisFunctionType, ResultType>::MaxwellIdentityOperator(
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label,
        int symmetry) :
    Base(domain, range, dualToRange, label, symmetry)
{
    typedef Fiber::MaxwellTestTrialIntegrandFunctor<
            BasisFunctionType, ResultType> IntegrandFunctor;
    typedef Fiber::DefaultTestTrialIntegral<IntegrandFunctor> Integral;
    m_integral.reset(new Integral((IntegrandFunctor())));
}

template <typename BasisFunctionType, typename ResultType>
MaxwellIdentityOperator<BasisFunctionType, ResultType>::~MaxwellIdentityOperator()
{
}

template <typename BasisFunctionType, typename ResultType>
const typename MaxwellIdentityOperator<BasisFunctionType, ResultType>::CollectionOfBasisTransformations&
MaxwellIdentityOperator<BasisFunctionType, ResultType>::testTransformations() const
{
    return this->dualToRange()->shapeFunctionValue();
}

template <typename BasisFunctionType, typename ResultType>
const typename MaxwellIdentityOperator<BasisFunctionType, ResultType>::CollectionOfBasisTransformations&
MaxwellIdentityOperator<BasisFunctionType, ResultType>::trialTransformations() const
{
    return this->domain()->shapeFunctionValue();
}

template <typename BasisFunctionType, typename ResultType>
const typename MaxwellIdentityOperator<BasisFunctionType, ResultType>::TestTrialIntegral&
MaxwellIdentityOperator<BasisFunctionType, ResultType>::integral() const
{
    return *m_integral;
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
maxwellIdentityOperator(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& domain,
                 const shared_ptr<const Space<BasisFunctionType> >& range,
                 const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
                 const std::string& label,
                 int symmetry)
{
    typedef MaxwellIdentityOperator<BasisFunctionType, ResultType> Id;
    return BoundaryOperator<BasisFunctionType, ResultType>(
                context,
                boost::make_shared<Id>(domain, range, dualToRange,
                                       label, symmetry));
}

#define INSTANTIATE_NONMEMBER_CONSTRUCTOR(BASIS, RESULT) \
    template BoundaryOperator<BASIS, RESULT> \
    maxwellIdentityOperator( \
        const shared_ptr<const Context<BASIS, RESULT> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const std::string&, \
        int)
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(MaxwellIdentityOperator);

} // namespace Bempp
