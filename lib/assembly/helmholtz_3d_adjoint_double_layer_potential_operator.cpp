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

#include "helmholtz_3d_adjoint_double_layer_potential_operator.hpp"
#include "helmholtz_3d_operator_base_imp.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/modified_helmholtz_3d_adjoint_double_layer_potential_kernel_functor.hpp"
#include "../fiber/scalar_function_value_functor.hpp"
#include "../fiber/simple_test_scalar_kernel_trial_integrand_functor.hpp"

#include "../fiber/standard_collection_of_kernels.hpp"
#include "../fiber/standard_collection_of_basis_transformations.hpp"
#include "../fiber/standard_test_kernel_trial_integral.hpp"

namespace Bempp
{

template <typename BasisFunctionType>
struct Helmholtz3dAdjointDoubleLayerPotentialOperatorImpl
{
    typedef Helmholtz3dAdjointDoubleLayerPotentialOperatorImpl<BasisFunctionType> This;
    typedef Helmholtz3dOperatorBase<This, BasisFunctionType> OperatorBase;
    typedef typename OperatorBase::CoordinateType CoordinateType;
    typedef typename OperatorBase::KernelType KernelType;
    typedef typename OperatorBase::ResultType ResultType;

    typedef Fiber::ModifiedHelmholtz3dAdjointDoubleLayerPotentialKernelFunctor<KernelType>
    KernelFunctor;
    typedef Fiber::ScalarFunctionValueFunctor<CoordinateType>
    TransformationFunctor;
    typedef Fiber::SimpleTestScalarKernelTrialIntegrandFunctor<
    BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

    explicit Helmholtz3dAdjointDoubleLayerPotentialOperatorImpl(KernelType waveNumber) :
        kernels(KernelFunctor(waveNumber / KernelType(0., 1.))),
        transformations(TransformationFunctor()),
        integral(IntegrandFunctor())
    {}

    Fiber::StandardCollectionOfKernels<KernelFunctor> kernels;
    Fiber::StandardCollectionOfBasisTransformations<TransformationFunctor>
    transformations;
    Fiber::StandardTestKernelTrialIntegral<IntegrandFunctor> integral;
};

template <typename BasisFunctionType>
Helmholtz3dAdjointDoubleLayerPotentialOperator<BasisFunctionType>::
Helmholtz3dAdjointDoubleLayerPotentialOperator(
        const Space<BasisFunctionType>& domain,
        const Space<BasisFunctionType>& range,
        const Space<BasisFunctionType>& dualToRange,
        KernelType waveNumber,
        const std::string& label) :
    Base(domain, range, dualToRange, waveNumber, label)
{
}

template <typename BasisFunctionType>
std::auto_ptr<LinearOperator<BasisFunctionType,
typename Helmholtz3dAdjointDoubleLayerPotentialOperator<BasisFunctionType>::ResultType> >
Helmholtz3dAdjointDoubleLayerPotentialOperator<BasisFunctionType>::
clone() const
{
    typedef LinearOperator<BasisFunctionType, ResultType> LinOp;
    typedef Helmholtz3dAdjointDoubleLayerPotentialOperator<BasisFunctionType> This;
    return std::auto_ptr<LinOp>(new This(*this));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dAdjointDoubleLayerPotentialOperator);

} // namespace Bempp
