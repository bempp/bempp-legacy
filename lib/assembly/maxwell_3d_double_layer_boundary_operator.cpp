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

#include "maxwell_3d_double_layer_boundary_operator.hpp"

#include "general_elementary_singular_integral_operator_imp.hpp"

#include "../common/boost_make_shared_fwd.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/modified_maxwell_3d_double_layer_potential_kernel_functor.hpp"
#include "../fiber/modified_maxwell_3d_double_layer_potential_kernel_interpolated_functor.hpp"
#include "../fiber/modified_maxwell_3d_double_layer_potential_integrand_functor.hpp"
#include "../fiber/hdiv_function_value_functor.hpp"

#include "../fiber/default_collection_of_kernels.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../fiber/default_test_kernel_trial_integral.hpp"

#include "../grid/max_distance.hpp"

namespace Bempp
{

template <typename BasisFunctionType>
BoundaryOperator<BasisFunctionType,
    typename ScalarTraits<BasisFunctionType>::ComplexType>
maxwell3dDoubleLayerBoundaryOperator(
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

    typedef Fiber::ModifiedMaxwell3dDoubleLayerPotentialKernelFunctor<KernelType>
        KernelFunctor;
    typedef Fiber::ModifiedMaxwell3dDoubleLayerPotentialKernelInterpolatedFunctor<KernelType>
        KernelInterpolatedFunctor;
    typedef Fiber::HdivFunctionValueFunctor<CoordinateType>
        TransformationFunctor;
    typedef Fiber::ModifiedMaxwell3dDoubleLayerPotentialIntegrandFunctor<
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
   maxwell3dDoubleLayerBoundaryOperator( \
       const shared_ptr<const Context<BASIS, ScalarTraits<BASIS>::ComplexType> >&, \
       const shared_ptr<const Space<BASIS> >&, \
       const shared_ptr<const Space<BASIS> >&, \
       const shared_ptr<const Space<BASIS> >&, \
       ScalarTraits<BASIS>::ComplexType, \
       const std::string&, int, bool, int)
FIBER_ITERATE_OVER_BASIS_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

} // namespace Bempp
