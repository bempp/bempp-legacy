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

#include "laplace_3d_single_layer_boundary_operator.hpp"
#include "laplace_3d_boundary_operator_base_imp.hpp"

#include "synthetic_scalar_integral_operator_builder.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/laplace_3d_single_layer_potential_kernel_functor.hpp"
#include "../fiber/scalar_function_value_functor.hpp"
#include "../fiber/simple_test_scalar_kernel_trial_integrand_functor.hpp"

#include "../fiber/default_collection_of_kernels.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../fiber/default_test_kernel_trial_integral.hpp"

#include "../common/boost_make_shared_fwd.hpp"

namespace Bempp
{

/** \cond PRIVATE */
template <typename BasisFunctionType, typename ResultType>
struct Laplace3dSingleLayerBoundaryOperatorImpl
{
    typedef Laplace3dSingleLayerBoundaryOperatorImpl<BasisFunctionType, ResultType>
    This;
    typedef Laplace3dBoundaryOperatorBase<This, BasisFunctionType, ResultType> BoundaryOperatorBase;
    typedef typename BoundaryOperatorBase::KernelType KernelType;
    typedef typename BoundaryOperatorBase::CoordinateType CoordinateType;

    typedef Fiber::Laplace3dSingleLayerPotentialKernelFunctor<KernelType>
    KernelFunctor;
    typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
    TransformationFunctor;
    typedef Fiber::SimpleTestScalarKernelTrialIntegrandFunctor<
    BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

    Laplace3dSingleLayerBoundaryOperatorImpl() :
        kernels(KernelFunctor()),
        transformations(TransformationFunctor()),
        integral(IntegrandFunctor())
    {}

    Fiber::DefaultCollectionOfKernels<KernelFunctor> kernels;
    Fiber::DefaultCollectionOfBasisTransformations<TransformationFunctor>
    transformations;
    Fiber::DefaultTestKernelTrialIntegral<IntegrandFunctor> integral;
};
/** \endcond */

template <typename BasisFunctionType, typename ResultType>
Laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType>::
Laplace3dSingleLayerBoundaryOperator(
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label,
        int symmetry) :
    Base(domain, range, dualToRange, label, symmetry)
{
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplace3dSingleLayerBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label,
        int symmetry)
{
   typedef Laplace3dSingleLayerBoundaryOperator<BasisFunctionType, ResultType> Op;
   return BoundaryOperator<BasisFunctionType, ResultType>(
               context, boost::make_shared<Op>(domain, range, dualToRange,
                                               label, symmetry));
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplace3dSyntheticSingleLayerBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const shared_ptr<const Space<BasisFunctionType> >& internalTrialSpace,
        const shared_ptr<const Space<BasisFunctionType> >& internalTestSpace,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY)
{
    BoundaryOperator<BasisFunctionType, ResultType> internalOp =
            laplace3dSingleLayerBoundaryOperator(
                context, internalTrialSpace, range /* or whatever */,
                internalTestSpace,
                "(" + label + ")_internal", symmetry);
    return makeSyntheticScalarIntegralOperator(
                internalOp, domain, range, dualToRange,
                internalTrialSpace, internalTestSpace,
                label, NO_SYMMETRY);
}

#define INSTANTIATE_NONMEMBER_CONSTRUCTOR(BASIS, RESULT) \
    template BoundaryOperator<BASIS, RESULT> \
    laplace3dSingleLayerBoundaryOperator( \
        const shared_ptr<const Context<BASIS, RESULT> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const std::string&, int); \
    template BoundaryOperator<BASIS, RESULT> \
    laplace3dSyntheticSingleLayerBoundaryOperator( \
        const shared_ptr<const Context<BASIS, RESULT> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const std::string&, int)
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

#define INSTANTIATE_BASE(BASIS, RESULT) \
    template class Laplace3dBoundaryOperatorBase< \
        Laplace3dSingleLayerBoundaryOperatorImpl<BASIS, RESULT>, BASIS, RESULT>
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_BASE);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dSingleLayerBoundaryOperator);

} // namespace Bempp
