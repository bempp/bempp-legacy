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
#include "ahmed_aux.hpp"
#include "../common/not_implemented_error.hpp"

#include <iostream>

namespace Bempp {

template <typename ValueType, typename GeometryTypeRows, typename GeometryTypeCols>
class DiscreteAcaScalarValuedLinearOperator :
        public DiscreteScalarValuedLinearOperator<ValueType>
{
public:
    DiscreteAcaScalarValuedLinearOperator(
            unsigned int rowCount, unsigned int columnCount,
            std::auto_ptr<bemblcluster<GeometryTypeRows, GeometryTypeCols> > blockCluster,
            mblock<ValueType>** blocks) :
        m_matrix(rowCount, columnCount, blockCluster, blocks)
    {}

    virtual void multiplyAddVector(ValueType multiplier,
                           const arma::Col<ValueType>& argument,
                           arma::Col<ValueType>& result)
    {
        // const_cast because Ahmed internally calls BLAS
        // functions, which don't respect const-correctness
        if (m_matrix.n != argument.n_rows)
            throw std::invalid_argument("DiscreteAcaScalarValuedLinearOperator::"
                                        "multiplyAddVector: "
                                        "incorrect vector length");
        result.set_size(argument.n_rows);
        m_matrix.amux(multiplier, const_cast<ValueType*>(argument.memptr()),
                      result.memptr());
    }

    virtual void dump() const
    {
        std::cout << asMatrix() << std::endl;
    }

    virtual arma::Mat<ValueType> asMatrix() const
    {
        const int rowCount = m_matrix.m;
        const int colCount = m_matrix.n;

        arma::Mat<ValueType> result(rowCount, colCount);
        result.fill(0.);
        arma::Col<ValueType> unit(colCount);
        unit.fill(0.);
        for (int col = 0; col < colCount; ++col)
        {
            if (col > 0)
                unit(col - 1) = 0.;
            unit(col) = 1.;
            m_matrix.amux(1., unit.memptr(), result.colptr(col));
        }
        return result;
    }

private:
    AhmedMatrix<ValueType, GeometryTypeRows, GeometryTypeCols> m_matrix;
};

} // namespace Bempp

#endif
