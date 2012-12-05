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

#ifndef bempp_general_elementary_singular_integral_operator_imp_hpp
#define bempp_general_elementary_singular_integral_operator_imp_hpp

#include "general_elementary_singular_integral_operator.hpp"

#include "../fiber/default_collection_of_kernels.hpp"
#include "../fiber/default_collection_of_basis_transformations.hpp"
#include "../fiber/default_test_kernel_trial_integral.hpp"

namespace Bempp
{

template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
template <typename KernelFunctor,
          typename TestTransformationsFunctor,
          typename TrialTransformationsFunctor,
          typename IntegrandFunctor>
GeneralElementarySingularIntegralOperator<
BasisFunctionType_, KernelType_, ResultType_>::
GeneralElementarySingularIntegralOperator(
        const shared_ptr<const Space<BasisFunctionType_> >& domain,
        const shared_ptr<const Space<BasisFunctionType_> >& range,
        const shared_ptr<const Space<BasisFunctionType_> >& dualToRange,
        const std::string& label,
        int symmetry,
        const KernelFunctor& kernelFunctor,
        const TestTransformationsFunctor& testTransformationsFunctor,
        const TrialTransformationsFunctor& trialTransformationsFunctor,
        const IntegrandFunctor& integrandFunctor) :
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
            integrandFunctor))
{
}

//template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
//const typename GeneralElementarySingularIntegralOperator<
//    BasisFunctionType_, KernelType_, ResultType_>::CollectionOfKernels&
//GeneralElementarySingularIntegralOperator<
//BasisFunctionType_, KernelType_, ResultType_>::kernels() const
//{
//    return *m_kernels;
//}

//template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
//const typename GeneralElementarySingularIntegralOperator<
//    BasisFunctionType_, KernelType_, ResultType_>::CollectionOfBasisTransformations&
//GeneralElementarySingularIntegralOperator<
//BasisFunctionType_, KernelType_, ResultType_>::testTransformations() const
//{
//    return *m_testTransformations;
//}

//template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
//const typename GeneralElementarySingularIntegralOperator<
//    BasisFunctionType_, KernelType_, ResultType_>::CollectionOfBasisTransformations&
//GeneralElementarySingularIntegralOperator<
//BasisFunctionType_, KernelType_, ResultType_>::trialTransformations() const
//{
//    return *m_trialTransformations;
//}

//template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
//const typename GeneralElementarySingularIntegralOperator<
//    BasisFunctionType_, KernelType_, ResultType_>::TestKernelTrialIntegral&
//GeneralElementarySingularIntegralOperator<
//BasisFunctionType_, KernelType_, ResultType_>::integral() const
//{
//    return *m_integral;
//}

//template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
//shared_ptr<const AbstractBoundaryOperatorId>
//GeneralElementarySingularIntegralOperator<
//BasisFunctionType_, KernelType_, ResultType_>::id() const
//{
//    return m_id;
//}

} // namespace Bempp

#endif
