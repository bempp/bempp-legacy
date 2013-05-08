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

#include "helmholtz_3d_hypersingular_boundary_operator.hpp"
#include "helmholtz_3d_boundary_operator_base_imp.hpp"

#include "general_hypersingular_integral_operator_imp.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/modified_helmholtz_3d_hypersingular_kernel_functor.hpp"
#include "../fiber/modified_helmholtz_3d_hypersingular_kernel_interpolated_functor.hpp"
#include "../fiber/modified_helmholtz_3d_hypersingular_off_diagonal_kernel_functor.hpp"
#include "../fiber/modified_helmholtz_3d_hypersingular_off_diagonal_interpolated_kernel_functor.hpp"
#include "../fiber/modified_helmholtz_3d_hypersingular_transformation_functor.hpp"
#include "../fiber/modified_helmholtz_3d_hypersingular_integrand_functor_2.hpp"
#include "../fiber/scalar_function_value_functor.hpp"
#include "../fiber/simple_test_scalar_kernel_trial_integrand_functor.hpp"

#include "../common/boost_make_shared_fwd.hpp"

namespace Bempp
{

template <typename BasisFunctionType>
BoundaryOperator<BasisFunctionType,
typename ScalarTraits<BasisFunctionType>::ComplexType>
helmholtz3dHypersingularBoundaryOperator(
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
    typedef typename ScalarTraits<BasisFunctionType>::ComplexType KernelType;
    typedef typename ScalarTraits<BasisFunctionType>::ComplexType ResultType;
    typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;

    typedef Fiber::ModifiedHelmholtz3dHypersingularKernelFunctor<KernelType>
            NoninterpolatedKernelFunctor;
    typedef Fiber::ModifiedHelmholtz3dHypersingularKernelInterpolatedFunctor<KernelType>
            InterpolatedKernelFunctor;
    typedef Fiber::ModifiedHelmholtz3dHypersingularTransformationFunctor<CoordinateType>
            TransformationFunctor;
    typedef Fiber::ModifiedHelmholtz3dHypersingularIntegrandFunctor2<
            BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

    typedef Fiber::ModifiedHelmholtz3dHypersingularOffDiagonalInterpolatedKernelFunctor<KernelType>
            OffDiagonalInterpolatedKernelFunctor;
    typedef Fiber::ModifiedHelmholtz3dHypersingularOffDiagonalKernelFunctor<KernelType>
            OffDiagonalNoninterpolatedKernelFunctor;
    typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
            OffDiagonalTransformationFunctor;
    typedef Fiber::SimpleTestScalarKernelTrialIntegrandFunctor<
            BasisFunctionType, KernelType, ResultType>
            OffDiagonalIntegrandFunctor;

    CoordinateType maxDistance_ =
            static_cast<CoordinateType>(1.1) *
            maxDistance(*domain->grid(), *dualToRange->grid());

    typedef GeneralHypersingularIntegralOperator<
            BasisFunctionType, KernelType, ResultType> Op;
    shared_ptr<Op> newOp;
    if (useInterpolation)
        newOp.reset(new Op(
                        domain, range, dualToRange, label, symmetry,
                        InterpolatedKernelFunctor(
                            waveNumber / KernelType(0., 1.),
                            maxDistance_,
                            interpPtsPerWavelength),
                        TransformationFunctor(),
                        TransformationFunctor(),
                        IntegrandFunctor(),
                        OffDiagonalInterpolatedKernelFunctor(
                            waveNumber / KernelType(0., 1.),
                            maxDistance_,
                            interpPtsPerWavelength),
                        OffDiagonalTransformationFunctor(),
                        OffDiagonalTransformationFunctor(),
                        OffDiagonalIntegrandFunctor()));
    else
        newOp.reset(new Op(
                        domain, range, dualToRange, label, symmetry,
                        NoninterpolatedKernelFunctor(
                            waveNumber / KernelType(0., 1.)),
                        TransformationFunctor(),
                        TransformationFunctor(),
                        IntegrandFunctor(),
                        OffDiagonalNoninterpolatedKernelFunctor(
                            waveNumber / KernelType(0., 1.)),
                        OffDiagonalTransformationFunctor(),
                        OffDiagonalTransformationFunctor(),
                        OffDiagonalIntegrandFunctor()));
    return BoundaryOperator<BasisFunctionType, ResultType>(context, newOp);
}

#define INSTANTIATE_NONMEMBER_CONSTRUCTOR(BASIS) \
   template BoundaryOperator<BASIS, ScalarTraits<BASIS>::ComplexType> \
   helmholtz3dHypersingularBoundaryOperator( \
       const shared_ptr<const Context<BASIS, ScalarTraits<BASIS>::ComplexType> >&, \
       const shared_ptr<const Space<BASIS> >&, \
       const shared_ptr<const Space<BASIS> >&, \
       const shared_ptr<const Space<BASIS> >&, \
       ScalarTraits<BASIS>::ComplexType, \
       const std::string&, int, bool, int)
FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

} // namespace Bempp
