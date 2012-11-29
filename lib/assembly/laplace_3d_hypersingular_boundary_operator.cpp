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

#include "laplace_3d_hypersingular_boundary_operator.hpp"
#include "laplace_3d_boundary_operator_base_imp.hpp"
#include "discrete_boundary_operator.hpp"

#include "context.hpp"
#include "identity_operator.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/laplace_3d_single_layer_potential_kernel_functor.hpp"
#include "../fiber/surface_curl_3d_functor.hpp"
#include "../fiber/simple_test_scalar_kernel_trial_integrand_functor.hpp"

#include "../fiber/default_collection_of_kernels.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../fiber/default_test_kernel_trial_integral.hpp"

#include "../space/unit_scalar_space.hpp"

namespace Bempp
{

/** \cond PRIVATE */
template <typename BasisFunctionType, typename ResultType>
struct Laplace3dHypersingularBoundaryOperatorImpl
{
    typedef typename Laplace3dBoundaryOperatorBase<BasisFunctionType, ResultType>::KernelType
    KernelType;
    typedef typename Laplace3dBoundaryOperatorBase<BasisFunctionType, ResultType>::CoordinateType
    CoordinateType;

    typedef Fiber::Laplace3dSingleLayerPotentialKernelFunctor<KernelType>
    KernelFunctor;
    typedef Fiber::SurfaceCurl3dFunctor<CoordinateType>
    TransformationFunctor;
    typedef Fiber::SimpleTestScalarKernelTrialIntegrandFunctor<
    BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

    Laplace3dHypersingularBoundaryOperatorImpl() :
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
Laplace3dHypersingularBoundaryOperator<BasisFunctionType, ResultType>::
Laplace3dHypersingularBoundaryOperator(
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
laplace3dHypersingularBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label,
        int symmetry)
{
    typedef Laplace3dHypersingularBoundaryOperator<BasisFunctionType, ResultType> Op;
    return BoundaryOperator<BasisFunctionType, ResultType>(
                context, boost::make_shared<Op>(domain, range, dualToRange,
                                                label, symmetry));
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplace3dModifiedHypersingularBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        ResultType alpha,
        const std::string& label,
        int symmetry)
{
    if (domain->grid() != range->grid())
        throw std::invalid_argument(
                "laplace3dModifiedHypersingularBoundaryOperator(): "
                "domain and range must be defined on the same grid");
    BoundaryOperator<BasisFunctionType, ResultType> hypOp =
            laplace3dHypersingularBoundaryOperator(
                context, domain, range, dualToRange, label, symmetry);

    AssemblyOptions idAssemblyOptions = context->assemblyOptions();
    idAssemblyOptions.enableSparseStorageOfMassMatrices(false);
    shared_ptr<Context<BasisFunctionType, ResultType> > idContext(
                new Context<BasisFunctionType, ResultType>(context->quadStrategy(),
                                                           idAssemblyOptions));
    shared_ptr<Space<BasisFunctionType> > unitSpace(
                new UnitScalarSpace<BasisFunctionType>(domain->grid()));
    BoundaryOperator<BasisFunctionType, ResultType> domainProjection =
            identityOperator<BasisFunctionType, ResultType>(
                idContext, domain, unitSpace, unitSpace, "DtU");
    BoundaryOperator<BasisFunctionType, ResultType> rangeProjection =
            identityOperator<BasisFunctionType, ResultType>(
                idContext, unitSpace, range, dualToRange, "UtR");

    BoundaryOperator<BasisFunctionType, ResultType> ii =
            identityOperator<BasisFunctionType, ResultType>(
                idContext, domain, range, dualToRange, "UtR");
    return hypOp + alpha * rangeProjection * domainProjection;
}

#define INSTANTIATE_NONMEMBER_CONSTRUCTOR(BASIS, RESULT) \
    template BoundaryOperator<BASIS, RESULT> \
    laplace3dHypersingularBoundaryOperator( \
    const shared_ptr<const Context<BASIS, RESULT> >&, \
    const shared_ptr<const Space<BASIS> >&, \
    const shared_ptr<const Space<BASIS> >&, \
    const shared_ptr<const Space<BASIS> >&, \
    const std::string&, int); \
    template BoundaryOperator<BASIS, RESULT> \
    laplace3dModifiedHypersingularBoundaryOperator( \
    const shared_ptr<const Context<BASIS, RESULT> >&, \
    const shared_ptr<const Space<BASIS> >&, \
    const shared_ptr<const Space<BASIS> >&, \
    const shared_ptr<const Space<BASIS> >&, \
    RESULT, const std::string&, int)
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

#define INSTANTIATE_BASE(BASIS, RESULT) \
    template class Laplace3dBoundaryOperatorBase< \
    Laplace3dHypersingularBoundaryOperatorImpl<BASIS, RESULT>, BASIS, RESULT>
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_BASE);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dHypersingularBoundaryOperator);

} // namespace Bempp
