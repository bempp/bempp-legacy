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

#include "laplace_3d_double_layer_potential_operator.hpp"
#include "laplace_3d_potential_operator_base_imp.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/laplace_3d_double_layer_potential_kernel_functor.hpp"
#include "../fiber/scalar_function_value_functor.hpp"
#include "../fiber/simple_scalar_kernel_trial_integrand_functor.hpp"

#include "../fiber/default_collection_of_kernels.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../fiber/default_kernel_trial_integral.hpp"

namespace Bempp
{

/** \cond PRIVATE */
template <typename BasisFunctionType, typename ResultType>
struct Laplace3dDoubleLayerPotentialOperatorImpl
{
    typedef Laplace3dDoubleLayerPotentialOperatorImpl<BasisFunctionType, ResultType>
    This;
    typedef Laplace3dPotentialOperatorBase<This, BasisFunctionType, ResultType> PotentialOperatorBase;
    typedef typename PotentialOperatorBase::KernelType KernelType;
    typedef typename PotentialOperatorBase::CoordinateType CoordinateType;

    typedef Fiber::Laplace3dDoubleLayerPotentialKernelFunctor<KernelType>
    KernelFunctor;
    typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
    TransformationFunctor;
    typedef Fiber::SimpleScalarKernelTrialIntegrandFunctor<
    BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

    Laplace3dDoubleLayerPotentialOperatorImpl() :
        kernels(KernelFunctor()),
        transformations(TransformationFunctor()),
        integral(IntegrandFunctor())
    {}

    Fiber::DefaultCollectionOfKernels<KernelFunctor> kernels;
    Fiber::DefaultCollectionOfBasisTransformations<TransformationFunctor>
    transformations;
    Fiber::DefaultKernelTrialIntegral<IntegrandFunctor> integral;
};
/** \endcond */

template <typename BasisFunctionType, typename ResultType>
Laplace3dDoubleLayerPotentialOperator<BasisFunctionType, ResultType>::
Laplace3dDoubleLayerPotentialOperator()
{
}

template <typename BasisFunctionType, typename ResultType>
Laplace3dDoubleLayerPotentialOperator<BasisFunctionType, ResultType>::
~Laplace3dDoubleLayerPotentialOperator()
{
}

#define INSTANTIATE_BASE_LAPLACE_DOUBLE_POTENTIAL(BASIS,RESULT)		     \
    template class Laplace3dPotentialOperatorBase< \
Laplace3dDoubleLayerPotentialOperatorImpl< BASIS , RESULT >, BASIS , RESULT >
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_BASE_LAPLACE_DOUBLE_POTENTIAL);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dDoubleLayerPotentialOperator);

} // namespace Bempp
