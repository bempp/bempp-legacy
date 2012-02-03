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

#ifndef bempp_concrete_quadrature_selector_hpp
#define bempp_concrete_quadrature_selector_hpp

#include "quadrature_selector.hpp"

#include "../common/types.hpp"

namespace Bempp
{

template <typename ValueType, typename FiberQuadratureSelector>
class ConcreteQuadratureSelector : public QuadratureSelector<ValueType>
{
private:
    FiberQuadratureSelector m_fiberQs;

public:
    explicit ConcreteQuadratureSelector(const FiberQuadratureSelector& fiberQs) :
        m_fiberQs(fiberQs) {
    }

    virtual void selectDoubleQuadratureRules(
            const std::vector<const GeometryAdapter*>& testGeometries,
            const std::vector<const GeometryAdapter*>& trialGeometries,
            Array2D<typename QuadratureRule>& quadRules) const {
        const int testGeometryCount = testGeometries.size();
        const int trialGeometryCount = trialGeometries.size();
        quadRules.set_size(testGeometryCount, trialGeometryCount);
        m_fiberQs.selectDoubleQuadratureRules(
                    testGeometryCount, testGeometries.begin(),
                    trialGeometryCount, trialGeometries.begin(),
                    quadRules.begin());
    }

    virtual void doubleQuadratureRulePointsAndWeights(
            typename QuadratureRule quadRule,
            arma::Mat<typename ctype>& testQuadPoints,
            arma::Mat<typename ctype>& trialQuadPoints,
            std::vector<ValueType>& quadWeights) const {
        const int pointCount =
                m_fiberQs.doubleQuadratureRulePointCount(quadRule);
        const int elementDim = m_fiberQs.elementDimension();
        testQuadPoints.set_size(elementDim, pointCount);
        trialQuadPoints.set_size(elementDim, pointCount);
        quadWeights.resize(pointCount);
        m_fiberQs.doubleQuadratureRulePointsAndWeights(
                    quadRule,
                    testQuadPoints.begin(),
                    trialQuadPoints.begin(),
                    quadWeights.begin());
    }
};

} // namespace Bempp


#endif
