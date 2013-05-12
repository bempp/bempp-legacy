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

#ifndef bempp_general_hypersingular_integral_operator_imp_hpp
#define bempp_general_hypersingular_integral_operator_imp_hpp

#include "general_hypersingular_integral_operator.hpp"

#include "../fiber/default_collection_of_kernels.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../fiber/default_test_kernel_trial_integral.hpp"

namespace Bempp
{

template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
template <typename KernelFunctor,
          typename TestTransformationsFunctor,
          typename TrialTransformationsFunctor,
          typename IntegrandFunctor,
          typename OffDiagonalKernelFunctor,
          typename OffDiagonalTestTransformationsFunctor,
          typename OffDiagonalTrialTransformationsFunctor,
          typename OffDiagonalIntegrandFunctor>
GeneralHypersingularIntegralOperator<
BasisFunctionType_, KernelType_, ResultType_>::
GeneralHypersingularIntegralOperator(
        const shared_ptr<const Space<BasisFunctionType_> >& domain,
        const shared_ptr<const Space<BasisFunctionType_> >& range,
        const shared_ptr<const Space<BasisFunctionType_> >& dualToRange,
        const std::string& label,
        int symmetry,
        const KernelFunctor& kernelFunctor,
        const TestTransformationsFunctor& testTransformationsFunctor,
        const TrialTransformationsFunctor& trialTransformationsFunctor,
        const IntegrandFunctor& integrandFunctor,
        const OffDiagonalKernelFunctor& offDiagonalKernelFunctor,
        const OffDiagonalTestTransformationsFunctor& offDiagonalTestTransformationsFunctor,
        const OffDiagonalTrialTransformationsFunctor& offDiagonalTrialTransformationsFunctor,
        const OffDiagonalIntegrandFunctor& offDiagonalIntegrandFunctor) :
    Base(domain, range, dualToRange, label, symmetry),
    m_kernels(
        new Fiber::DefaultCollectionOfKernels<KernelFunctor>(kernelFunctor)),
    m_testTransformations(
        new Fiber::DefaultCollectionOfBasisTransformations<TestTransformationsFunctor>(
            testTransformationsFunctor)),
    m_trialTransformations(
        new Fiber::DefaultCollectionOfBasisTransformations<TrialTransformationsFunctor>(
            trialTransformationsFunctor)),
    m_integral(
        new Fiber::DefaultTestKernelTrialIntegral<IntegrandFunctor>(
            integrandFunctor)),
    m_offDiagonalKernels(
        new Fiber::DefaultCollectionOfKernels<OffDiagonalKernelFunctor>(
            offDiagonalKernelFunctor)),
    m_offDiagonalTestTransformations(
        new Fiber::DefaultCollectionOfBasisTransformations<OffDiagonalTestTransformationsFunctor>(
            offDiagonalTestTransformationsFunctor)),
    m_offDiagonalTrialTransformations(
        new Fiber::DefaultCollectionOfBasisTransformations<OffDiagonalTrialTransformationsFunctor>(
            offDiagonalTrialTransformationsFunctor)),
    m_offDiagonalIntegral(
        new Fiber::DefaultTestKernelTrialIntegral<OffDiagonalIntegrandFunctor>(
            offDiagonalIntegrandFunctor))
{
}

} // namespace Bempp

#endif
