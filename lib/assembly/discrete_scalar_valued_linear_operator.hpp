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

#ifndef bempp_discrete_scalar_valued_linear_operator_hpp
#define bempp_discrete_scalar_valued_linear_operator_hpp

#include <armadillo>

namespace Bempp {

template <typename ValueType>
class DiscreteScalarValuedLinearOperator
{
public:
    virtual ~DiscreteScalarValuedLinearOperator() {}

    /** \brief Matrix-vector multiplication.

      Sets \p result to

      result := result + multiplier * L * argument,

      where L is the linear operator represented by this object. */
    virtual void multiplyAddVector(ValueType multiplier,
                           const arma::Col<ValueType>& argument,
                           arma::Col<ValueType>& result) = 0;

    /** \brief Write a textual representation of the operator to standard output. */
    virtual void dump() const = 0;

    /** \brief Matrix representation of the operator.

    This method need not be supported by all subclasses. */
    virtual arma::Mat<ValueType> asMatrix() const = 0;

    /** \brief Calculate a block of this operator's matrix representation and
      add it to \p block.

      \param[in] rows  Vector of row indices.
      \param[in] cols  Vector of col indices.
      \param[in,out] block
        On entry, matrix of size (rows.size(), cols.size()).
        On exit, each element (i, j) of this matrix will be augmented by the
        element (rows[i], cols[j]) of the matrix representation of this operator.

      Row and column indices may be unsorted.

      This method need not be supported by all subclasses. */
    virtual void addBlock(const std::vector<int>& rows,
                          const std::vector<int>& cols,
                          arma::Mat<ValueType>& block) const = 0;
};

} // namespace Bempp

#endif
