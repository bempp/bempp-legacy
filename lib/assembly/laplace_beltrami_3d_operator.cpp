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

#include "laplace_beltrami_3d_operator.hpp"

#include "boundary_operator.hpp"

#include "../common/boost_make_shared_fwd.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../fiber/default_test_trial_integral_imp.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/simple_test_trial_integrand_functor.hpp"
#include "../fiber/surface_grad_3d_functor.hpp"

#include <boost/type_traits/is_complex.hpp>

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
LaplaceBeltrami3dOperator<BasisFunctionType, ResultType>::LaplaceBeltrami3dOperator(
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label,
        int symmetry) :
    Base(domain, range, dualToRange, label, symmetry)
{
    if (domain->codomainDimension() != 1)
        throw std::invalid_argument(
                "LaplaceBeltrami3dOperator::LaplaceBeltrami3dOperator(): "
                "domain must consist of scalar-valued functions");
    if (range->codomainDimension() != 1)
        throw std::invalid_argument(
                "LaplaceBeltrami3dOperator::LaplaceBeltrami3dOperator(): "
                "range must consist of scalar-valued functions");
    if (dualToRange->codomainDimension() != 1)
        throw std::invalid_argument(
                "LaplaceBeltrami3dOperator::LaplaceBeltrami3dOperator(): "
                "space dual to range must consist of scalar-valued functions");

    typedef Fiber::SimpleTestTrialIntegrandFunctor<
            BasisFunctionType, ResultType> IntegrandFunctor;
    typedef Fiber::DefaultTestTrialIntegral<IntegrandFunctor> Integral;
    typedef Fiber::SurfaceGrad3dFunctor<CoordinateType> TransformationFunctor;
    typedef Fiber::DefaultCollectionOfBasisTransformations<TransformationFunctor> Transformation;
    m_surfaceGrad.reset(new Transformation((TransformationFunctor())));
    m_integral.reset(new Integral((IntegrandFunctor())));
}

template <typename BasisFunctionType, typename ResultType>
LaplaceBeltrami3dOperator<BasisFunctionType, ResultType>::~LaplaceBeltrami3dOperator()
{
}

template <typename BasisFunctionType, typename ResultType>
const typename LaplaceBeltrami3dOperator<BasisFunctionType, ResultType>::CollectionOfBasisTransformations&
LaplaceBeltrami3dOperator<BasisFunctionType, ResultType>::testTransformations() const
{
    return *m_surfaceGrad;
}

template <typename BasisFunctionType, typename ResultType>
const typename LaplaceBeltrami3dOperator<BasisFunctionType, ResultType>::CollectionOfBasisTransformations&
LaplaceBeltrami3dOperator<BasisFunctionType, ResultType>::trialTransformations() const
{
    return *m_surfaceGrad;
}

template <typename BasisFunctionType, typename ResultType>
const typename LaplaceBeltrami3dOperator<BasisFunctionType, ResultType>::TestTrialIntegral&
LaplaceBeltrami3dOperator<BasisFunctionType, ResultType>::integral() const
{
    return *m_integral;
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplaceBeltrami3dOperator(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& domain,
                 const shared_ptr<const Space<BasisFunctionType> >& range,
                 const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
                 const std::string& label,
                 int symmetry)
{
    typedef LaplaceBeltrami3dOperator<BasisFunctionType, ResultType> Id;
    return BoundaryOperator<BasisFunctionType, ResultType>(
                context,
                boost::make_shared<Id>(domain, range, dualToRange,
                                       label, symmetry));
}

#define INSTANTIATE_NONMEMBER_CONSTRUCTOR(BASIS, RESULT) \
    template BoundaryOperator<BASIS, RESULT> \
    laplaceBeltrami3dOperator( \
        const shared_ptr<const Context<BASIS, RESULT> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const std::string&, \
        int)
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(LaplaceBeltrami3dOperator);

} // namespace Bempp
