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

#include "laplace_3d_hypersingular_operator.hpp"
#include "laplace_3d_operator_base_imp.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../fiber/laplace_3d_single_layer_potential_kernel_functor.hpp"
#include "../fiber/surface_curl_3d_functor.hpp"
#include "../fiber/simple_test_scalar_kernel_trial_integrand_functor.hpp"

#include "../fiber/standard_collection_of_kernels.hpp"
#include "../fiber/standard_collection_of_basis_transformations.hpp"
#include "../fiber/standard_test_kernel_trial_integral.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
struct Laplace3dHypersingularOperatorImpl
{
    typedef typename Laplace3dOperatorBase<BasisFunctionType, ResultType>::KernelType
    KernelType;
    typedef typename Laplace3dOperatorBase<BasisFunctionType, ResultType>::CoordinateType
    CoordinateType;

    typedef Fiber::Laplace3dSingleLayerPotentialKernelFunctor<KernelType>
    KernelFunctor;
    typedef Fiber::SurfaceCurl3dFunctor<CoordinateType>
    TransformationFunctor;
    typedef Fiber::SimpleTestScalarKernelTrialIntegrandFunctor<
    BasisFunctionType, KernelType, ResultType> IntegrandFunctor;

    Laplace3dHypersingularOperatorImpl() :
        kernels(KernelFunctor()),
        transformations(TransformationFunctor()),
        integral(IntegrandFunctor())
    {}

    Fiber::StandardCollectionOfKernels<KernelFunctor> kernels;
    Fiber::StandardCollectionOfBasisTransformations<TransformationFunctor>
    transformations;
    Fiber::StandardTestKernelTrialIntegral<IntegrandFunctor> integral;
};

template <typename BasisFunctionType, typename ResultType>
Laplace3dHypersingularOperator<BasisFunctionType, ResultType>::
Laplace3dHypersingularOperator(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace) :
    Base(testSpace, trialSpace)
{
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Laplace3dHypersingularOperator);

} // namespace Bempp
