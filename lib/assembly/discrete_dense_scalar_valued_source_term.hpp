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

#ifndef bempp_discrete_dense_scalar_valued_source_term_hpp
#define bempp_discrete_dense_scalar_valued_source_term_hpp

#include "discrete_scalar_valued_source_term.hpp"

#include <armadillo>
#include <iostream>

namespace Bempp {

template <typename ValueType>
class DiscreteDenseScalarValuedSourceTerm :
        public DiscreteScalarValuedSourceTerm<ValueType>
{
public:
    DiscreteDenseScalarValuedSourceTerm(const arma::Col<ValueType>& vec) :
        m_vec(vec) {}

    virtual void dump() const {
        std::cout << m_vec << std::endl;
    }

    virtual arma::Col<ValueType> asVector() const {
        return m_vec;
    }

private:
    arma::Col<ValueType> m_vec;
};

} // namespace Bempp

#endif
