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
            boost::ptr_vector<TermType>& terms,
            const std::vector<ValueType>& multipliers) :
        m_multipliers(multipliers)
    {
        m_terms.transfer(m_terms.end(), terms);
        if (m_terms.size() != m_multipliers.size())
            throw std::invalid_argument(
                    "DiscreteScalarValuedLinearOperatorSuperposition::"
                    "DiscreteScalarValuedLinearOperatorSuperposition(): "
                    "incompatible argument lengths");
    }

    virtual void multiplyAddVector(ValueType multiplier,
                           const arma::Col<ValueType>& argument,
                           arma::Col<ValueType>& result)
    {
        for (int i = 0; i < m_terms.size(); ++i)
            m_terms[i].multiplyAddVector(m_multipliers[i], argument, result);
    }

    virtual void dump() const
    {
        throw NotImplementedError(
                    "DiscreteScalarValuedLinearOperatorSuperposition::dump(): "
                    "not implemented yet");
    }

private:
    boost::ptr_vector<TermType> m_terms;
    std::vector<ValueType> m_multipliers;
};

} // namespace Bempp

#endif
