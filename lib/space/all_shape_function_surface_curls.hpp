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

#ifndef bempp_all_shape_function_surface_curls_hpp
#define bempp_all_shape_function_surface_curls_hpp

#include "function_family.hpp"
#include "space.hpp"

#include <armadillo>
#include <vector>

namespace Bempp {

template <typename ValueType>
class AllShapeFunctionSurfaceCurls : public FunctionFamily<ValueType>
{
public:
    AllShapeFunctionSurfaceCurls(const Space<ValueType>& space,
                                 ElementVariant elementVariant) :
        m_space(space), m_elementVariant(elementVariant)
    {}

    virtual int domainDimension() const {
        return m_space.domainDimension();
    }

    virtual int codomainDimension() const {
        return -1; // FIXME: what is the correct value? world dim. or grid dim.?
    }

    virtual bool needsJacobianInverseTransposed() const {
        return true;
    }

    virtual int size() const {
        return m_space.basisFunctionCount(m_elementVariant);
    }

    // coordCount might be different from domainDimension() if
    // e.g. barycentric coordinates are used
    // localQuadPoints[coordCount,pointCount] (column-first)
    virtual void setEvaluationPoints(const arma::Mat<ctype>& evaluationPoints) {
        m_evaluationPoints = evaluationPoints;
        m_space.evaluateBasisFunctions(m_elementVariant, m_evaluationPoints,
                                       m_basisFunctions);
        m_basisFunctionDerivatives.resize(domainDimension());
        for (int i = 0; i < domainDimension(); ++i)
            m_space.evaluateBasisFunctionDerivative(
                        m_elementVariant, m_evaluationPoints, i,
                        m_basisFunctionDerivatives[i]);
    }

    virtual int evaluationPointCount() const {
        return m_evaluationPoints.n_cols;
    }

    // result[codomainDimension,functionCount(),pointCount]
    virtual void evaluate(const arma::Cube<ctype>& jacobianInverseTransposed,
                          arma::Cube<ValueType>& result) const {
        m_space.evaluateShapeFunctionSurfaceCurlsInternal(
                    m_basisFunctions, m_basisFunctionDerivatives,
                    jInvT, res);
    }

private:
    const Space<ValueType>& m_space;
    int m_elementVariant;
    arma::Mat<ctype> m_localQuadPoints;
    arma::Cube<ValueType> m_basisFunctions;
    std::vector<arma::Cube<ValueType> > m_basisFunctionDerivatives;
};

// namespace Bempp

#endif
