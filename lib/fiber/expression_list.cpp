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

#include "expression_list.hpp"
#include "explicit_instantiation.hpp"

namespace Fiber
{

template <typename ResultType>
size_t ExpressionList<ResultType>::domainDimension() const
{
    if (m_terms.empty())
        return 0;
    else
        return m_terms[0]->domainDimension();
}

template <typename ResultType>
size_t ExpressionList<ResultType>::codomainDimension() const
{
    if (m_terms.empty())
        return 0;
    else
        return m_terms[0]->codomainDimension();
}

template <typename ResultType>
void ExpressionList<ResultType>::addDependencies(
        size_t& basisDeps, size_t& geomDeps) const
{
    for (size_t i = 0; i < m_terms.size(); ++i)
        m_terms[i]->addDependencies(basisDeps, geomDeps);
}

template <typename ResultType>
void ExpressionList<ResultType>::addTerm(
        const Expression<CoordinateType>& expression,
        ResultType weight)
{
    if (!m_terms.empty())
        if (expression.domainDimension() != domainDimension() ||
                expression.codomainDimension() != codomainDimension())
            throw std::invalid_argument("ExpressionList::addTerm(): "
                                        "domain or codomain dimension incompatible "
                                        "with those of the other terms");
    m_terms.push_back(&expression);
    m_weights.push_back(weight);
}

template <typename ResultType>
bool ExpressionList<ResultType>::isTrivial() const
{
    return (m_terms.size() == 1 &&
            m_weights[0] == static_cast<ResultType>(1.));
}

template <typename ResultType>
void ExpressionList<ResultType>::clear()
{
    m_terms.clear();
    m_weights.clear();
}

template <typename ResultType>
void ExpressionList<ResultType>::evaluate(
        const BasisData<CoordinateType>& basisData,
        const GeometricalData<CoordinateType>& geomData,
        std::vector<arma::Cube<CoordinateType> >& result) const
{
    const size_t termCount = m_terms.size();
    result.resize(termCount);
    for (size_t i = 0; i < termCount; ++i)
        m_terms[i]->evaluate(basisData, geomData, result[i]);
}

template <typename ResultType>
void ExpressionList<ResultType>::evaluate(
        const BasisData<ComplexType>& basisData,
        const GeometricalData<CoordinateType>& geomData,
        std::vector<arma::Cube<ComplexType> >& result) const
{
    const size_t termCount = m_terms.size();
    result.resize(termCount);
    for (size_t i = 0; i < termCount; ++i)
        m_terms[i]->evaluate(basisData, geomData, result[i]);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(ExpressionList);

} // namespace Fiber

