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
            mblock<ValueType>** blocks,
            const std::vector<unsigned int>& permutedRowIndices,
            const std::vector<unsigned int>& originalColumnIndices) :
        m_matrix(rowCount, columnCount, blockCluster, blocks),
        m_permutedRowIndices(permutedRowIndices),
        m_originalColumnIndices(originalColumnIndices)
    {}

    /** \brief Matrix-vector multiplication.

      Sets \p result to

      result := result + multiplier * L * argument,

      where L is the linear operator represented by this object.

      \internal The vectors \p argument and \p result are in "standard" ordering.
      This method performs automatically the necessary conversions to and from
      the permuted ordering used by Ahmed. */
    virtual void multiplyAddVector(ValueType multiplier,
                                   const arma::Col<ValueType>& argument,
                                   arma::Col<ValueType>& result)
    {
        arma::Col<ValueType> permutedArgument;
        permuteArgument(argument, permutedArgument);

        // const_cast because Ahmed internally calls BLAS
        // functions, which don't respect const-correctness
        if (m_matrix.n != argument.n_rows)
            throw std::invalid_argument("DiscreteAcaScalarValuedLinearOperator::"
                                        "multiplyAddVector: "
                                        "incorrect vector length");
        arma::Col<ValueType> permutedResult;
        permutedResult.set_size(argument.n_rows);
        m_matrix.amux(multiplier,
                      const_cast<ValueType*>(permutedArgument.memptr()),
                      permutedResult.memptr());
        unpermuteResult(permutedResult, result);
    }

    virtual void dump() const
    {
        std::cout << asMatrix() << std::endl;
    }

    virtual arma::Mat<ValueType> asMatrix() const
    {
        const int rowCount = m_matrix.m;
        const int colCount = m_matrix.n;

        arma::Mat<ValueType> permutedOutput(rowCount, colCount);
        permutedOutput.fill(0.);
        arma::Col<ValueType> unit(colCount);
        unit.fill(0.);

        for (int col = 0; col < colCount; ++col)
        {
            if (col > 0)
                unit(col - 1) = 0.;
            unit(col) = 1.;
            m_matrix.amux(1., unit.memptr(), permutedOutput.colptr(col));
        }

        arma::Mat<ValueType> output(rowCount, colCount);
        for (int col = 0; col < colCount; ++col)
            for (int row = 0; row < rowCount; ++row)
                output(row, m_originalColumnIndices[col]) =
                        permutedOutput(m_permutedRowIndices[row], col);

        return output;
    }

    /** \brief Convert an argument of the operator from original to permuted ordering. */
    void permuteArgument(const arma::Col<ValueType>& argument,
                         arma::Col<ValueType>& permutedArgument) const
    {
        const int dim = argument.n_elem;
        permutedArgument.set_size(dim);
        for (int i = 0; i < dim; ++i)
            permutedArgument(i) = argument(m_originalColumnIndices[i]);
    }

    /** \brief Convert a result of the operator's action from permuted to original ordering. */
    void unpermuteResult(const arma::Col<ValueType>& permutedResult,
                         arma::Col<ValueType>& result) const
    {
        const int dim = permutedResult.n_elem;
        result.set_size(dim);
        for (int i = 0; i < dim; ++i)
            result(i) = permutedResult(m_permutedRowIndices[i]);
    }

private:
    AhmedMatrix<ValueType, GeometryTypeRows, GeometryTypeCols> m_matrix;
    std::vector<unsigned int> m_permutedRowIndices;
    std::vector<unsigned int> m_originalColumnIndices;
};

} // namespace Bempp

#endif
