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

#ifndef bempp_elementary_singular_integral_operator_hpp
#define bempp_elementary_singular_integral_operator_hpp

#include "../common/common.hpp"

#include "elementary_integral_operator.hpp"

namespace Bempp
{

template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
class ElementarySingularIntegralOperator :
        public ElementaryIntegralOperator<BasisFunctionType_, KernelType_, ResultType_>
{
    typedef ElementaryIntegralOperator<BasisFunctionType_, KernelType_, ResultType_> Base;
public:
    typedef typename Base::BasisFunctionType BasisFunctionType;
    typedef typename Base::KernelType KernelType;
    typedef typename Base::ResultType ResultType;
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::CollectionOfBasisTransformations
    CollectionOfBasisTransformations;
    typedef typename Base::CollectionOfKernels CollectionOfKernels;
    typedef typename Base::TestKernelTrialIntegral TestKernelTrialIntegral;

    ElementarySingularIntegralOperator(
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
            const std::string& label = "") :
        Base(domain, range, dualToRange, label) {
    }

    virtual bool isRegular() const {
        return false;
    }
};

} // namespace Bempp

#endif
