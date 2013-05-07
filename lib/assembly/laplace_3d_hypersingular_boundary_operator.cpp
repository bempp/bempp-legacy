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
#include "general_elementary_local_operator_imp.hpp"
#include "general_hypersingular_integral_operator_imp.hpp"
#include "identity_operator.hpp"
#include "laplace_3d_single_layer_boundary_operator.hpp"
#include "synthetic_integral_operator.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/laplace_3d_single_layer_potential_kernel_functor.hpp"
#include "../fiber/laplace_3d_hypersingular_off_diagonal_kernel_functor.hpp"
#include "../fiber/surface_curl_3d_functor.hpp"
#include "../fiber/scalar_function_value_functor.hpp"
#include "../fiber/simple_test_scalar_kernel_trial_integrand_functor.hpp"
#include "../fiber/single_component_test_trial_integrand_functor.hpp"

#include "../fiber/default_collection_of_kernels.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../fiber/default_test_kernel_trial_integral.hpp"

#include "../space/unit_scalar_space.hpp"

#include <boost/type_traits/is_complex.hpp>

namespace Bempp
{

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
    typedef typename ScalarTraits<BasisFunctionType>::RealType KernelType;
    typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

    typedef Fiber::Laplace3dSingleLayerPotentialKernelFunctor<KernelType>
    KernelFunctor;
    typedef Fiber::SurfaceCurl3dFunctor<CoordinateType>
    TransformationFunctor;
    typedef Fiber::SimpleTestScalarKernelTrialIntegrandFunctor<
    BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

    typedef Fiber::Laplace3dHypersingularOffDiagonalKernelFunctor<KernelType>
    OffDiagonalKernelFunctor;
    typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
    OffDiagonalTransformationFunctor;
    typedef Fiber::SimpleTestScalarKernelTrialIntegrandFunctor<
    BasisFunctionType, KernelType, ResultType> OffDiagonalIntegrandFunctor;

    typedef GeneralHypersingularIntegralOperator<
            BasisFunctionType, KernelType, ResultType> Op;
    shared_ptr<Op> newOp(new Op(
                             domain, range, dualToRange, label, symmetry,
                             KernelFunctor(),
                             TransformationFunctor(),
                             TransformationFunctor(),
                             IntegrandFunctor(),
                             OffDiagonalKernelFunctor(),
                             OffDiagonalTransformationFunctor(),
                             OffDiagonalTransformationFunctor(),
                             OffDiagonalIntegrandFunctor()));
    return BoundaryOperator<BasisFunctionType, ResultType>(context, newOp);
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplace3dSyntheticHypersingularBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const shared_ptr<const Space<BasisFunctionType> >& internalTrialSpace,
        const shared_ptr<const Space<BasisFunctionType> >& internalTestSpace,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY)
{
    // Note: we don't really need to care about ranges and duals to domains of
    // the internal operator. The only range space that matters is that of the
    // leftmost operator in the product.

    const char xyz[] = "xyz";
    const size_t dimWorld = 3;

    BoundaryOperator<BasisFunctionType, ResultType> slp =
            laplace3dSingleLayerBoundaryOperator(
                context, internalTrialSpace, internalTestSpace /* or whatever */,
                internalTestSpace,
                "(" + label + ")_internal", symmetry);

    typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

    typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
            ValueFunctor;
    typedef Fiber::SurfaceCurl3dFunctor<CoordinateType>
            CurlFunctor;
    typedef Fiber::SingleComponentTestTrialIntegrandFunctor<
            BasisFunctionType, ResultType> IntegrandFunctor;

    typedef GeneralElementaryLocalOperator<BasisFunctionType, ResultType> LocalOp;

    std::vector<BoundaryOperator<BasisFunctionType, ResultType> > trialCurlComponents;
    std::vector<BoundaryOperator<BasisFunctionType, ResultType> > testCurlComponents;
        testCurlComponents.resize(3);
        for (size_t i = 0; i < dimWorld; ++i)
            testCurlComponents[i] = BoundaryOperator<BasisFunctionType, ResultType>(
                        context, boost::make_shared<LocalOp>(
                            internalTestSpace, range, dualToRange,
                            ("(" + label + ")_test_curl_") + xyz[i], NO_SYMMETRY,
                            CurlFunctor(),
                            ValueFunctor(),
                            IntegrandFunctor(i, 0)));
    size_t overallSymmetry = 0; // symmetry of the decomposition
    if (domain == dualToRange && internalTrialSpace == internalTestSpace)
        overallSymmetry = HERMITIAN |
                (boost::is_complex<BasisFunctionType>() ? 0 : SYMMETRIC);
    else {
        trialCurlComponents.resize(3);
        for (size_t i = 0; i < dimWorld; ++i)
            trialCurlComponents[i] = BoundaryOperator<BasisFunctionType, ResultType>(
                        context, boost::make_shared<LocalOp>(
                            domain, internalTrialSpace /* or whatever */, internalTrialSpace,
                            ("(" + label + ")_trial_curl_") + xyz[i], NO_SYMMETRY,
                            ValueFunctor(),
                            CurlFunctor(),
                            IntegrandFunctor(0, i)));
    }

    typedef SyntheticIntegralOperator<BasisFunctionType, ResultType> SyntheticOp;

    return BoundaryOperator<BasisFunctionType, ResultType>(
                context, boost::make_shared<SyntheticOp>(
                    testCurlComponents, slp, trialCurlComponents, label,
                    overallSymmetry));
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
    laplace3dSyntheticHypersingularBoundaryOperator( \
    const shared_ptr<const Context<BASIS, RESULT> >&, \
    const shared_ptr<const Space<BASIS> >&, \
    const shared_ptr<const Space<BASIS> >&, \
    const shared_ptr<const Space<BASIS> >&, \
    const shared_ptr<const Space<BASIS> >&, \
    const shared_ptr<const Space<BASIS> >&, \
    const std::string&, int);
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

} // namespace Bempp
