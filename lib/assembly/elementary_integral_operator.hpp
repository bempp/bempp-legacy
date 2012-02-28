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
#include "../fiber/kernel.hpp"

#include "../common/types.hpp"

#include <memory>
#include <vector>

namespace Fiber
{

template <typename ValueType> class Expression;

} // namespace Fiber

namespace Bempp
{

class Geometry;

template <typename ValueType>
class ElementaryIntegralOperator : public ElementaryLinearOperator<ValueType>
{
public:
    typedef typename ElementaryLinearOperator<ValueType>::IntegrationManager
    IntegrationManager;
    typedef typename ElementaryLinearOperator<ValueType>::IntegrationManagerFactory
    IntegrationManagerFactory;

    virtual int trialComponentCount() const {
        return kernel().domainDimension();
    }

    virtual int testComponentCount() const {
        return kernel().codomainDimension();
    }

    virtual bool isRegular() const = 0;

private:
    virtual std::auto_ptr<IntegrationManager > makeIntegrationManager(
            const IntegrationManagerFactory& factory,
            const GeometryFactory& geometryFactory,
            const arma::Mat<ValueType>& vertices,
            const arma::Mat<int>& elementCorners,
            const arma::Mat<char>& auxData) const;

    virtual const Fiber::Kernel<ValueType>& kernel() const = 0;
    virtual const Fiber::Expression<ValueType>& testExpression() const = 0;
    virtual const Fiber::Expression<ValueType>& trialExpression() const = 0;

    virtual void evaluateLocalWeakForms(
            CallVariant evalVariant,
            const std::vector<const EntityPointer<0>*>& elementsA,
            const EntityPointer<0>& elementB,
            LocalDofIndex localDofIndexB,
            const Space<ValueType>& spaceA,
            const Space<ValueType>& spaceB,
            IntegrationManager& intMgr,
            std::vector<arma::Mat<ValueType> >& result) const;

    virtual void evaluateLocalWeakForms(
            const std::vector<const EntityPointer<0>*>& testElements,
            const std::vector<const EntityPointer<0>*>& trialElements,
            const Space<ValueType>& testSpace,
            const Space<ValueType>& trialSpace,
            IntegrationManager& intMgr,
            Fiber::Array2D<arma::Mat<ValueType> >& result) const;
};

} // namespace Bempp

#endif
