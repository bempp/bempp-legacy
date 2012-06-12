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

#include "helmholtz_3d_hypersingular_operator.hpp"
#include "helmholtz_3d_operator_base_imp.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/modified_helmholtz_3d_single_layer_potential_kernel_functor.hpp"
#include "../fiber/modified_helmholtz_3d_hypersingular_transformation_functor.hpp"
#include "../fiber/modified_helmholtz_3d_hypersingular_integrand_functor.hpp"

#include "../fiber/standard_collection_of_kernels.hpp"
#include "../fiber/standard_collection_of_basis_transformations.hpp"
#include "../fiber/standard_test_kernel_trial_integral.hpp"

//    m_expressionList.addTerm(m_surfaceCurl);
//    m_expressionList.addTerm(m_valueTimesNormal, KernelType(0., 1.) * waveNumber);

namespace Bempp
{

template <typename BasisFunctionType>
struct Helmholtz3dHypersingularOperatorImpl
{
    typedef Helmholtz3dHypersingularOperatorImpl<BasisFunctionType> This;
    typedef Helmholtz3dOperatorBase<This, BasisFunctionType> OperatorBase;
    typedef typename OperatorBase::CoordinateType CoordinateType;
    typedef typename OperatorBase::KernelType KernelType;
    typedef typename OperatorBase::ResultType ResultType;

    typedef Fiber::ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<KernelType>
    KernelFunctor;
    typedef Fiber::ModifiedHelmholtz3dHypersingularTransformationFunctor<CoordinateType>
    TransformationFunctor;
    typedef Fiber::ModifiedHelmholtz3dHypersingularIntegrandFunctor<
    BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

    explicit Helmholtz3dHypersingularOperatorImpl(KernelType waveNumber) :
        kernels(KernelFunctor(waveNumber / KernelType(0., 1.))),
        transformations(TransformationFunctor()),
        integral(IntegrandFunctor(waveNumber / KernelType(0., 1.)))
    {}

    Fiber::StandardCollectionOfKernels<KernelFunctor> kernels;
    Fiber::StandardCollectionOfBasisTransformations<TransformationFunctor>
    transformations;
    Fiber::StandardTestKernelTrialIntegral<IntegrandFunctor> integral;
};

template <typename BasisFunctionType>
Helmholtz3dHypersingularOperator<BasisFunctionType>::
Helmholtz3dHypersingularOperator(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        KernelType waveNumber) :
    Base(testSpace, trialSpace, waveNumber)
{
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(Helmholtz3dHypersingularOperator);

} // namespace Bempp
