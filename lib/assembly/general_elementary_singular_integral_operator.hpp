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

#ifndef bempp_general_elementary_singular_integral_operator_hpp
#define bempp_general_elementary_singular_integral_operator_hpp

#include "elementary_singular_integral_operator.hpp"

namespace Bempp
{

template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
class GeneralElementarySingularIntegralOperator :
        public ElementarySingularIntegralOperator<
            BasisFunctionType_, KernelType_, ResultType_>
{
    typedef ElementarySingularIntegralOperator<
    BasisFunctionType_, KernelType_, ResultType_> Base;
public:
    /** \brief Type of the values of the basis functions into which functions
     *  acted upon by the operator are expanded. */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \brief Type of the values of kernel functions. */
    typedef typename Base::KernelType KernelType;
    /** \copydoc ElementaryIntegralOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc ElementaryIntegralOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc ElementaryIntegralOperator::CollectionOfBasisTransformations */
    typedef typename Base::CollectionOfBasisTransformations
    CollectionOfBasisTransformations;
    /** \copydoc ElementaryIntegralOperator::CollectionOfKernels */
    typedef typename Base::CollectionOfKernels CollectionOfKernels;
    /** \copydoc ElementaryIntegralOperator::TestKernelTrialIntegral */
    typedef typename Base::TestKernelTrialIntegral TestKernelTrialIntegral;

    template <typename KernelFunctor,
              typename TestTransformationsFunctor,
              typename TrialTransformationsFunctor,
              typename IntegrandFunctor>
    GeneralElementarySingularIntegralOperator(
            const shared_ptr<const Space<BasisFunctionType_> >& domain,
            const shared_ptr<const Space<BasisFunctionType_> >& range,
            const shared_ptr<const Space<BasisFunctionType_> >& dualToRange,
            const std::string& label,
            int symmetry,
            const KernelFunctor& kernelFunctor,
            const TestTransformationsFunctor& testTransformationsFunctor,
            const TrialTransformationsFunctor& trialTransformationsFunctor,
            const IntegrandFunctor& integrandFunctor);

    virtual const CollectionOfKernels& kernels() const
    { return *m_kernels; }
    virtual const CollectionOfBasisTransformations& testTransformations() const
    { return *m_testTransformations; }
    virtual const CollectionOfBasisTransformations& trialTransformations() const
    { return *m_trialTransformations; }
    virtual const TestKernelTrialIntegral& integral() const
    { return *m_integral; }

private:
    /** \cond PRIVATE */
    shared_ptr<CollectionOfKernels> m_kernels;
    shared_ptr<CollectionOfBasisTransformations> m_testTransformations;
    shared_ptr<CollectionOfBasisTransformations> m_trialTransformations;
    shared_ptr<TestKernelTrialIntegral> m_integral;
    /** \endcond */
};

} // namespace Bempp

#endif
