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

#ifndef fiber_expression_list_hpp
#define fiber_expression_list_hpp

#include "../common/common.hpp"

#include "expression.hpp"
#include "scalar_traits.hpp"

#include "../common/armadillo_fwd.hpp"
#include <vector>

namespace Fiber
{

template <typename ResultType>
class ExpressionList
{
public:
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
    typedef typename ScalarTraits<ResultType>::ComplexType ComplexType;

    int domainDimension() const;
    int codomainDimension() const;

    void addDependencies(size_t& basisDeps, size_t& geomDeps) const;
    void addTerm(const Expression<CoordinateType>& expression,
                 ResultType weight = ResultType(1.));
    size_t termCount() const { return m_terms.size(); }
    /** \brief Return true if the list consists of only one term with unit weight. */
    bool isTrivial() const;
    void clear();

    void evaluate(const BasisData<CoordinateType>& basisData,
                       const GeometricalData<CoordinateType>& geomData,
                       std::vector<arma::Cube<CoordinateType> >& result) const;
    void evaluate(const BasisData<ComplexType>& basisData,
                       const GeometricalData<CoordinateType>& geomData,
                       std::vector<arma::Cube<ComplexType> >& result) const;

    const Expression<CoordinateType>& term(size_t i) const { return *m_terms[i]; }
    ResultType weight(size_t i) const { return m_weights[i]; }

private:
    std::vector<const Expression<CoordinateType>*> m_terms;
    std::vector<ResultType> m_weights;
};

} // namespace Fiber

#endif // fiber_expression_list_hpp
