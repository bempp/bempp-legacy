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

#ifndef bempp_discrete_scalar_valued_linear_operator_superposition_hpp
#define bempp_discrete_scalar_valued_linear_operator_superposition_hpp

#include "discrete_scalar_valued_linear_operator.hpp"
#include "../common/not_implemented_error.hpp"

#include <boost/ptr_container/ptr_vector.hpp>

namespace Bempp {

template <typename ValueType>
class DiscreteScalarValuedLinearOperatorSuperposition :
        public DiscreteScalarValuedLinearOperator<ValueType>
{
public:
    typedef DiscreteScalarValuedLinearOperator<ValueType> TermType;

    /* acquires ownership of these operators */
    DiscreteScalarValuedLinearOperatorSuperposition(
            boost::ptr_vector<TermType>& terms)
    {
        m_terms.transfer(m_terms.end(), terms);
    }

    virtual void multiplyAddVector(ValueType multiplier,
                           const arma::Col<ValueType>& argument,
                           arma::Col<ValueType>& result)
    {
        for (int i = 0; i < m_terms.size(); ++i)
            m_terms[i].multiplyAddVector(multiplier, argument, result);
    }

    virtual void dump() const
    {
        throw NotImplementedError(
                    "DiscreteScalarValuedLinearOperatorSuperposition::dump(): "
                    "not implemented yet");
    }

    virtual arma::Mat<ValueType> asMatrix() const
    {
        arma::Mat<ValueType> result;
        if (!m_terms.empty()) {
            result = m_terms[0].asMatrix();
            for (int i = 1; i < m_terms.size(); ++i)
                result += m_terms[i].asMatrix();
        }
        return result;
    }

    virtual void addBlock(const std::vector<int>& rows,
                          const std::vector<int>& cols,
                          arma::Mat<ValueType>& block) const
    {
        // TODO: implement this
        throw std::runtime_error("DiscreteScalarValuedLinearOperatorSuperposition::"
                                 "addBlock(): Not implemented yet");
        for (int i = 0; i < m_terms.size(); ++i)
                m_terms[i].addBlock(rows, cols, block);
    }

private:
    boost::ptr_vector<TermType> m_terms;
};

} // namespace Bempp

#endif
