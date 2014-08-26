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

namespace Bempp {

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
template <typename KernelFunctor, typename TestTransformationsFunctor,
          typename TrialTransformationsFunctor, typename IntegrandFunctor>
GeneralElementarySingularIntegralOperator<BasisFunctionType_, KernelType_,
                                          ResultType_>::
    GeneralElementarySingularIntegralOperator(
        const shared_ptr<const Space<BasisFunctionType_>> &domain,
        const shared_ptr<const Space<BasisFunctionType_>> &range,
        const shared_ptr<const Space<BasisFunctionType_>> &dualToRange,
        const std::string &label, int symmetry,
        const KernelFunctor &kernelFunctor,
        const TestTransformationsFunctor &testTransformationsFunctor,
        const TrialTransformationsFunctor &trialTransformationsFunctor,
        const IntegrandFunctor &integrandFunctor)
    : Base(domain, range, dualToRange, label, symmetry),
      m_kernels(
          new Fiber::DefaultCollectionOfKernels<KernelFunctor>(kernelFunctor)),
      m_testTransformations(
          new Fiber::DefaultCollectionOfShapesetTransformations<
              TestTransformationsFunctor>(testTransformationsFunctor)),
      m_trialTransformations(
          new Fiber::DefaultCollectionOfShapesetTransformations<
              TrialTransformationsFunctor>(trialTransformationsFunctor)),
      m_integral(new Fiber::DefaultTestKernelTrialIntegral<IntegrandFunctor>(
          integrandFunctor)) {}

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
template <typename KernelFunctor, typename TestTransformationsFunctor,
          typename TrialTransformationsFunctor>
GeneralElementarySingularIntegralOperator<BasisFunctionType_, KernelType_,
                                          ResultType_>::
    GeneralElementarySingularIntegralOperator(
        const shared_ptr<const Space<BasisFunctionType_>> &domain,
        const shared_ptr<const Space<BasisFunctionType_>> &range,
        const shared_ptr<const Space<BasisFunctionType_>> &dualToRange,
        const std::string &label, int symmetry,
        const KernelFunctor &kernelFunctor,
        const TestTransformationsFunctor &testTransformationsFunctor,
        const TrialTransformationsFunctor &trialTransformationsFunctor,
        const shared_ptr<Fiber::TestKernelTrialIntegral<
            BasisFunctionType_, KernelType_, ResultType_>> &integral)
    : Base(domain, range, dualToRange, label, symmetry),
      m_kernels(
          new Fiber::DefaultCollectionOfKernels<KernelFunctor>(kernelFunctor)),
      m_testTransformations(
          new Fiber::DefaultCollectionOfShapesetTransformations<
              TestTransformationsFunctor>(testTransformationsFunctor)),
      m_trialTransformations(
          new Fiber::DefaultCollectionOfShapesetTransformations<
              TrialTransformationsFunctor>(trialTransformationsFunctor)),
      m_integral(integral) {}

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
GeneralElementarySingularIntegralOperator<BasisFunctionType_, KernelType_,
                                          ResultType_>::
    GeneralElementarySingularIntegralOperator(
        const shared_ptr<const Space<BasisFunctionType_>> &domain,
        const shared_ptr<const Space<BasisFunctionType_>> &range,
        const shared_ptr<const Space<BasisFunctionType_>> &dualToRange,
        const std::string &label, int symmetry,
        const shared_ptr<Fiber::CollectionOfKernels<KernelType_>> &kernels,
        const shared_ptr<Fiber::CollectionOfShapesetTransformations<
            CoordinateType>> &testTransformations,
        const shared_ptr<Fiber::CollectionOfShapesetTransformations<
            CoordinateType>> &trialTransformations,
        const shared_ptr<Fiber::TestKernelTrialIntegral<
            BasisFunctionType_, KernelType_, ResultType_>> &integral)
    : Base(domain, range, dualToRange, label, symmetry), m_kernels(kernels),
      m_testTransformations(testTransformations),
      m_trialTransformations(trialTransformations), m_integral(integral) {}

} // namespace Bempp

#endif
