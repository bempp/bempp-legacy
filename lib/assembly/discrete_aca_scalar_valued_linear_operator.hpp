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

#ifndef bempp_discrete_aca_scalar_valued_linear_operator_hpp
#define bempp_discrete_aca_scalar_valued_linear_operator_hpp

#include "discrete_scalar_valued_linear_operator.hpp"
#include "ahmed_matrix.hpp"

namespace Bempp {

template <typename ValueType, typename GeometryTypeRows, typename GeometryTypeCols>
class DiscreteAcaScalarValuedLinearOperator :
        public DiscreteScalarValuedLinearOperator<ValueType>
{
public:
    DiscreteAcaScalarValuedLinearOperator(
            unsigned int rowCount, unsigned int columnCount,
            auto_ptr<bemblcluster<GeometryTypeRows, GeometryTypeCols> > blockCluster,
            mblock<ValueType>** blocks) :
        m_matrix(rowCount, columnCount, blockCluster, blocks)
    {}

    virtual void multiplyAddVector(ValueType multiplier,
                           const arma::Col<ValueType>& argument,
                           arma::Col<ValueType>& result)
    {
        // TODO: ensure that the arguments have correct dimensions.
        m_matrix.amux(multiplier, argument.memptr(), result.memptr());
    }

private:
    AhmedMatrix m_matrix<ValueType, GeometryTypeRows, GeometryTypeCols>;
};

} // namespace Bempp

#endif
