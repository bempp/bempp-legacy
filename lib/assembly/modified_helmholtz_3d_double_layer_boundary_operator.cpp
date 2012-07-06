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

#include "modified_helmholtz_3d_double_layer_boundary_operator.hpp"
#include "modified_helmholtz_3d_boundary_operator_base_imp.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/modified_helmholtz_3d_double_layer_potential_kernel_functor.hpp"
#include "../fiber/scalar_function_value_functor.hpp"
#include "../fiber/simple_test_scalar_kernel_trial_integrand_functor.hpp"

#include "../fiber/default_collection_of_kernels.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../fiber/default_test_kernel_trial_integral.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename KernelType_, typename ResultType>
struct ModifiedHelmholtz3dDoubleLayerBoundaryOperatorImpl
{
    typedef KernelType_ KernelType;
    typedef ModifiedHelmholtz3dDoubleLayerBoundaryOperatorImpl<
    BasisFunctionType, KernelType, ResultType> This;
    typedef ModifiedHelmholtz3dBoundaryOperatorBase<This,
    BasisFunctionType, KernelType, ResultType> BoundaryOperatorBase;
    typedef typename BoundaryOperatorBase::CoordinateType CoordinateType;

    typedef Fiber::ModifiedHelmholtz3dDoubleLayerPotentialKernelFunctor<KernelType>
    KernelFunctor;
    typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
    TransformationFunctor;
    typedef Fiber::SimpleTestScalarKernelTrialIntegrandFunctor<
    BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

    explicit ModifiedHelmholtz3dDoubleLayerBoundaryOperatorImpl(KernelType waveNumber) :
        kernels(KernelFunctor(waveNumber)),
        transformations(TransformationFunctor()),
        integral(IntegrandFunctor())
    {}

    Fiber::DefaultCollectionOfKernels<KernelFunctor> kernels;
    Fiber::DefaultCollectionOfBasisTransformations<TransformationFunctor>
    transformations;
    Fiber::DefaultTestKernelTrialIntegral<IntegrandFunctor> integral;
};

template <typename BasisFunctionType, typename KernelType, typename ResultType>
ModifiedHelmholtz3dDoubleLayerBoundaryOperator<BasisFunctionType, KernelType, ResultType>::
ModifiedHelmholtz3dDoubleLayerBoundaryOperator(
        const Space<BasisFunctionType>& domain,
        const Space<BasisFunctionType>& range,
        const Space<BasisFunctionType>& dualToRange,
        KernelType waveNumber,
        const std::string& label) :
    Base(domain, range, dualToRange, waveNumber, label)
{
}

#define INSTANTIATE_BASE(BASIS, KERNEL, RESULT) \
    template class ModifiedHelmholtz3dBoundaryOperatorBase< \
    ModifiedHelmholtz3dDoubleLayerBoundaryOperatorImpl<BASIS, KERNEL, RESULT>, \
    BASIS, KERNEL, RESULT>
FIBER_ITERATE_OVER_BASIS_KERNEL_AND_RESULT_TYPES(INSTANTIATE_BASE);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(
        ModifiedHelmholtz3dDoubleLayerBoundaryOperator);

} // namespace Bempp
