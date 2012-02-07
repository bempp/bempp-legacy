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

#ifndef bempp_elementary_integral_operator_hpp
#define bempp_elementary_integral_operator_hpp

#include "elementary_linear_operator.hpp"
#include "kernel.hpp"

#include "../common/types.hpp"

#include <memory>

namespace Bempp {

template <typename ValueType> class FunctionFamily;

template <typename ValueType>
class ElementaryIntegralOperator : public ElementaryLinearOperator<ValueType>
{
public:
    virtual int trialComponentCount() const {
        return kernel().domainDimension();
    }

    virtual int testComponentCount() const {
        return kernel().codomainDimension();
    }

    virtual bool isRegular() const = 0;

private:
    virtual const Kernel<ValueType>& kernel() const = 0;
    virtual std::auto_ptr<FunctionFamily<ValueType> > testFunctionFamily(
            const Space<ValueType>& testSpace,
            ElementVariant elementVariant) const = 0;
    virtual std::auto_ptr<FunctionFamily<ValueType> > trialFunctionFamily(
            const Space<ValueType>& trialSpace,
            ElementVariant elementVariant) const = 0;

    virtual void evaluateLocalWeakForms(
            const std::vector<const EntityPointer<0>*>& testElements,
            const std::vector<const EntityPointer<0>*>& trialElements,
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            const QuadratureSelector<ValueType>& quadSelector,
            const AssemblyOptions& options,
            Array2D<arma::Mat<ValueType> >& result) const;
    void evaluateLocalWeakFormsWithOpenCl(
            const std::vector<const EntityPointer<0>*>& testElements,
            const std::vector<const EntityPointer<0>*>& trialElements,
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            const QuadratureSelector<ValueType>& quadSelector,
            const AssemblyOptions& options,
            Array2D<arma::Mat<ValueType> >& result) const;
    void evaluateLocalWeakFormsWithoutOpenCl(
            const std::vector<const EntityPointer<0>*>& testElements,
            const std::vector<const EntityPointer<0>*>& trialElements,
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            const QuadratureSelector<ValueType>& quadSelector,
            const AssemblyOptions& options,
            Array2D<arma::Mat<ValueType> >& result) const;
};

} // namespace Bempp

#endif
