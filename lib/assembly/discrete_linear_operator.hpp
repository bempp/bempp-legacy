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

#include "config_trilinos.hpp"

#ifndef bempp_discrete_linear_operator_hpp
#define bempp_discrete_linear_operator_hpp

#include <armadillo>

#ifdef WITH_TRILINOS
#include <Thyra_LinearOpDefaultBase_decl.hpp>
#endif

namespace Bempp {

enum TranspositionMode
{
#ifdef WITH_TRILINOS
    NO_TRANSPOSE = Thyra::NOTRANS,
    CONJUGATE = Thyra::CONJ,
    TRANSPOSE = Thyra::TRANS,
    CONJUGATE_TRANSPOSE = Thyra::CONJTRANS
#else
    NO_TRANSPOSE,
    CONJUGATE,
    TRANSPOSE,
    CONJUGATE_TRANSPOSE
#endif
};

template <typename ValueType>
class DiscreteLinearOperator
#ifdef WITH_TRILINOS
        : public Thyra::LinearOpDefaultBase<ValueType>
#endif
{
public:
    virtual ~DiscreteLinearOperator() {}

#ifdef WITH_TRILINOS
    // import the apply() member function from base class
    // before we overload it
    using Thyra::LinearOpDefaultBase<ValueType>::apply;
#endif

    /** \brief Matrix-vector multiplication.

      Sets

      \p y_inout := \p alpha * L * \p x_in + \p beta * \p y_inout,

      where L is the linear operator represented by this object.

      This overload is always supported, even if compiling without Trilinos. */
    void apply(const TranspositionMode trans,
               const arma::Col<ValueType>& x_in,
               arma::Col<ValueType>& y_inout,
               const ValueType alpha,
               const ValueType beta) const {
        applyBuiltInImpl(trans, x_in, y_inout, alpha, beta);
    }

    /** \brief Write a textual representation of the operator to standard output. */
    virtual void dump() const = 0;

    /** \brief Matrix representation of the operator.

    This method need not be supported by all subclasses. */
    virtual arma::Mat<ValueType> asMatrix() const = 0;

    /** \brief Number of rows of the operator. */
    virtual unsigned int rowCount() const = 0;

    /** \brief Number of columns of the operator. */
    virtual unsigned int columnCount() const = 0;

    /** \brief Perform the operation \p block += alpha*A[\p rows,\p cols], where
               \p alpha is a scalar and A[\rows,\p cols] is a subblock of the linear
               operator A.

      \param[in] rows  Vector of row indices.
      \param[in] cols  Vector of col indices.
      \param[in] alpha scalar multiplier for the block
      \param[in,out] block
        On entry, matrix of size (rows.size(), cols.size()).
        On exit, each element (i, j) of this matrix will be augmented by the
        element (rows[i], cols[j]) multiplied with alpha of the matrix representation of this operator.

      Row and column indices may be unsorted.

      This method need not be supported by all subclasses. */
    virtual void addBlock(const std::vector<int>& rows,
                          const std::vector<int>& cols,
                          const ValueType alpha,
                          arma::Mat<ValueType>& block) const = 0;

private:
    virtual void applyBuiltInImpl(const TranspositionMode trans,
                                  const arma::Col<ValueType>& x_in,
                                  arma::Col<ValueType>& y_inout,
                                  const ValueType alpha,
                                  const ValueType beta) const = 0;
};

} // namespace Bempp

#endif
