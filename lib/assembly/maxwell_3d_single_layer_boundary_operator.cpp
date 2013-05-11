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

#include "maxwell_3d_single_layer_boundary_operator.hpp"

#include "context.hpp"
#include "general_elementary_local_operator_imp.hpp"
#include "general_elementary_singular_integral_operator_imp.hpp"
#include "general_hypersingular_integral_operator_imp.hpp"
#include "helmholtz_3d_single_layer_boundary_operator.hpp"
#include "synthetic_integral_operator.hpp"

#include "../common/boost_make_shared_fwd.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/modified_maxwell_3d_single_layer_boundary_operator_kernel_functor.hpp"
#include "../fiber/modified_maxwell_3d_single_layer_boundary_operator_kernel_interpolated_functor.hpp"
#include "../fiber/modified_maxwell_3d_single_layer_operators_transformation_functor.hpp"
#include "../fiber/modified_maxwell_3d_single_layer_boundary_operator_integrand_functor.hpp"
#include "../fiber/hdiv_function_value_functor.hpp"
#include "../fiber/scalar_function_value_functor.hpp"
#include "../fiber/simple_test_trial_integrand_functor.hpp"
#include "../fiber/single_component_test_trial_integrand_functor.hpp"
#include "../fiber/surface_div_3d_functor.hpp"

#include "../fiber/default_collection_of_kernels.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../fiber/default_test_kernel_trial_integral.hpp"

#include "../grid/max_distance.hpp"

#include <boost/type_traits/is_complex.hpp>

namespace Bempp
{


template <typename BasisFunctionType>
BoundaryOperator<BasisFunctionType,
    typename ScalarTraits<BasisFunctionType>::ComplexType>
maxwell3dSyntheticSingleLayerBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType,
        typename ScalarTraits<BasisFunctionType>::ComplexType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        typename ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const std::string& label,
        int internalSymmetry,
        bool useInterpolation,
        int interpPtsPerWavelength)
{
    typedef typename ScalarTraits<BasisFunctionType>::ComplexType KernelType;
    typedef typename ScalarTraits<BasisFunctionType>::ComplexType ResultType;
    typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

    if (!domain || !range || !dualToRange)
        throw std::invalid_argument(
            "maxwell3dSyntheticSingleLayerBoundaryOperator(): "
            "domain, range and dualToRange must not be null");
    AssemblyOptions internalAssemblyOptions = context->assemblyOptions();
    AcaOptions internalAcaOptions = internalAssemblyOptions.acaOptions();
    internalAcaOptions.globalAssemblyBeforeCompression = true;
    internalAcaOptions.mode = AcaOptions::GLOBAL_ASSEMBLY;
    internalAssemblyOptions.switchToAcaMode(internalAcaOptions);
    typedef Context<BasisFunctionType, ResultType> Ctx;
    shared_ptr<Ctx> internalContext(new Ctx(
            context->quadStrategy(), internalAssemblyOptions));
    shared_ptr<const Space<BasisFunctionType> > internalTrialSpace = 
        domain->discontinuousSpace(domain);
    shared_ptr<const Space<BasisFunctionType> > internalTestSpace = 
        dualToRange->discontinuousSpace(dualToRange);

    typedef GeneralElementaryLocalOperator<BasisFunctionType, ResultType> LocalOp;
    typedef SyntheticIntegralOperator<BasisFunctionType, ResultType> SyntheticOp;

    // Note: we don't really need to care about ranges and duals to domains of
    // the internal operator. The only range space that matters is that of the
    // leftmost operator in the product.

    const char xyz[] = "xyz";
    const size_t dimWorld = 3;
    const ResultType kappa = waveNumber / ResultType(0, 1);
    const ResultType invKappa = static_cast<CoordinateType>(1.) / kappa;

    BoundaryOperator<BasisFunctionType, ResultType> slp =
            helmholtz3dSingleLayerBoundaryOperator<BasisFunctionType>(
                internalContext, internalTrialSpace, internalTestSpace /* or whatever */,
                internalTestSpace,
                waveNumber, "(" + label + ")_internal", internalSymmetry,
                useInterpolation, interpPtsPerWavelength);

    typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
            ScalarValueFunctor;
    typedef Fiber::SurfaceDiv3dFunctor<CoordinateType>
            DivFunctor;
    typedef Fiber::HdivFunctionValueFunctor<CoordinateType>
            VectorValueFunctor;
    typedef Fiber::SingleComponentTestTrialIntegrandFunctor<
            BasisFunctionType, ResultType> Term0IntegrandFunctor;
    typedef Fiber::SimpleTestTrialIntegrandFunctor<
            BasisFunctionType, ResultType> Term1IntegrandFunctor;

    // symmetry of the decomposition
    size_t syntheseSymmetry = 0;
    if (domain == dualToRange && internalTrialSpace == internalTestSpace)
        syntheseSymmetry = HERMITIAN |
                (boost::is_complex<BasisFunctionType>() ? 0 : SYMMETRIC);

    std::vector<BoundaryOperator<BasisFunctionType, ResultType> > testLocalOps;
    std::vector<BoundaryOperator<BasisFunctionType, ResultType> > trialLocalOps;
    testLocalOps.resize(dimWorld);
    for (size_t i = 0; i < dimWorld; ++i)
        testLocalOps[i] =
                BoundaryOperator<BasisFunctionType, ResultType>(
                    internalContext, boost::make_shared<LocalOp>(
                        internalTestSpace, range, dualToRange,
                        ("(" + label + ")_test_") + xyz[i], NO_SYMMETRY,
                        VectorValueFunctor(),
                        ScalarValueFunctor(),
                        Term0IntegrandFunctor(i, 0)));

    if (!syntheseSymmetry) {
        trialLocalOps.resize(dimWorld);
        for (size_t i = 0; i < dimWorld; ++i)
            trialLocalOps[i] =
                    BoundaryOperator<BasisFunctionType, ResultType>(
                        internalContext, boost::make_shared<LocalOp>(
                            domain, internalTrialSpace /* or whatever */,
                            internalTrialSpace,
                            ("(" + label + ")_trial_") + xyz[i], NO_SYMMETRY,
                            ScalarValueFunctor(),
                            VectorValueFunctor(),
                            Term0IntegrandFunctor(0, i)));
    }
    // It might be more prudent to distinguish between the symmetry of the total
    // operator and the symmetry of the decomposition
    BoundaryOperator<BasisFunctionType, ResultType> term0(
                context, boost::make_shared<SyntheticOp>(
                    testLocalOps, kappa * slp, trialLocalOps,
                    "(" + label + ")_term_1", syntheseSymmetry));

    testLocalOps.resize(1);
    testLocalOps[0] =
            BoundaryOperator<BasisFunctionType, ResultType>(
                internalContext, boost::make_shared<LocalOp>(
                    internalTestSpace, range, dualToRange,
                    ("(" + label + ")_test_div"), NO_SYMMETRY,
                    DivFunctor(),
                    ScalarValueFunctor(),
                    Term1IntegrandFunctor()));
    if (!syntheseSymmetry) {
        trialLocalOps.resize(1);
        trialLocalOps[0] =
                BoundaryOperator<BasisFunctionType, ResultType>(
                    internalContext, boost::make_shared<LocalOp>(
                        domain, internalTrialSpace /* or whatever */,
                        internalTrialSpace,
                        ("(" + label + ")_trial_div"), NO_SYMMETRY,
                        ScalarValueFunctor(),
                        DivFunctor(),
                        Term1IntegrandFunctor()));
    } else
        trialLocalOps.clear();
    BoundaryOperator<BasisFunctionType, ResultType> term1(
                context, boost::make_shared<SyntheticOp>(
                    testLocalOps, invKappa * slp, trialLocalOps,
                    "(" + label + ")_term_2", syntheseSymmetry));

    return term0 + term1;
}

template <typename BasisFunctionType>
BoundaryOperator<BasisFunctionType,
    typename ScalarTraits<BasisFunctionType>::ComplexType>
maxwell3dSingleLayerBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType,
        typename ScalarTraits<BasisFunctionType>::ComplexType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        typename ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
        const std::string& label,
        int symmetry,
        bool useInterpolation,
        int interpPtsPerWavelength)
{
    const AssemblyOptions& assemblyOptions = context->assemblyOptions();
    if (assemblyOptions.assemblyMode() == AssemblyOptions::ACA &&
        (!assemblyOptions.acaOptions().globalAssemblyBeforeCompression ||
         assemblyOptions.acaOptions().mode == AcaOptions::LOCAL_ASSEMBLY))
        return maxwell3dSyntheticSingleLayerBoundaryOperator(
            context, domain, range, dualToRange, waveNumber, label, 
            symmetry, useInterpolation, interpPtsPerWavelength);

    typedef typename ScalarTraits<BasisFunctionType>::ComplexType KernelType;
    typedef typename ScalarTraits<BasisFunctionType>::ComplexType ResultType;
    typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

    typedef Fiber::ModifiedMaxwell3dSingleLayerBoundaryOperatorKernelFunctor<KernelType>
        KernelFunctor;
    typedef Fiber::ModifiedMaxwell3dSingleLayerBoundaryOperatorKernelInterpolatedFunctor<KernelType>
        KernelInterpolatedFunctor;
    typedef Fiber::ModifiedMaxwell3dSingleLayerOperatorsTransformationFunctor<CoordinateType>
        TransformationFunctor;
    typedef Fiber::ModifiedMaxwell3dSingleLayerBoundaryOperatorIntegrandFunctor<
        BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

    typedef GeneralElementarySingularIntegralOperator<
            BasisFunctionType, KernelType, ResultType> Op;
    if (useInterpolation)
        return BoundaryOperator<BasisFunctionType, ResultType>(
                context, boost::make_shared<Op>(
                    domain, range, dualToRange, label, symmetry,
                    KernelInterpolatedFunctor(waveNumber / KernelType(0., 1.),
                                              1.1 * maxDistance(
                                                  *domain->grid(),
                                                  *dualToRange->grid()),
                                              interpPtsPerWavelength),
                    TransformationFunctor(),
                    TransformationFunctor(),
                    IntegrandFunctor()));
    else
        return BoundaryOperator<BasisFunctionType, ResultType>(
                context, boost::make_shared<Op>(
                    domain, range, dualToRange, label, symmetry,
                    KernelFunctor(waveNumber / KernelType(0., 1.)),
                    TransformationFunctor(),
                    TransformationFunctor(),
                    IntegrandFunctor()));
}

#define INSTANTIATE_NONMEMBER_CONSTRUCTOR(BASIS) \
    template BoundaryOperator<BASIS, ScalarTraits<BASIS>::ComplexType> \
    maxwell3dSingleLayerBoundaryOperator( \
        const shared_ptr<const Context<BASIS, ScalarTraits<BASIS>::ComplexType> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        ScalarTraits<BASIS>::ComplexType, \
        const std::string&, int, bool, int)
FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

} // namespace Bempp
