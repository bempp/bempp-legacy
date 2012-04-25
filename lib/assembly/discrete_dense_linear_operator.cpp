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

#include "discrete_dense_linear_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <iostream>
#include <stdexcept>

#ifdef WITH_TRILINOS
#include <Thyra_DetachedSpmdVectorView.hpp>
#include <Thyra_SpmdVectorSpaceDefaultBase.hpp>
#endif

namespace Bempp
{

template <typename ValueType>
DiscreteDenseLinearOperator<ValueType>::
DiscreteDenseLinearOperator(const arma::Mat<ValueType>& mat) :
    m_mat(mat)
#ifdef WITH_TRILINOS
  , m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(mat.n_cols)),
    m_rangeSpace(Thyra::defaultSpmdVectorSpace<ValueType>(mat.n_rows))
#endif
{
}

template <typename ValueType>
void DiscreteDenseLinearOperator<ValueType>::dump() const
{
    std::cout << m_mat << std::endl;
}

template <typename ValueType>
arma::Mat<ValueType> DiscreteDenseLinearOperator<ValueType>::asMatrix() const
{
    return m_mat;
}

template <typename ValueType>
unsigned int DiscreteDenseLinearOperator<ValueType>::rowCount() const
{
    return m_mat.n_rows;
}

template <typename ValueType>
unsigned int DiscreteDenseLinearOperator<ValueType>::columnCount() const
{
    return m_mat.n_cols;
}

template <typename ValueType>
void DiscreteDenseLinearOperator<ValueType>::addBlock(
        const std::vector<int>& rows,
        const std::vector<int>& cols,
        const ValueType alpha,
        arma::Mat<ValueType>& block) const
{
    if (block.n_rows != rows.size() || block.n_cols != cols.size())
        throw std::invalid_argument(
                "DiscreteDenseLinearOperator::addBlock(): "
                "incorrect block size");
    for (int col = 0; col < cols.size(); ++col)
        for (int row = 0; row < rows.size(); ++row)
            block(row, col) += alpha*m_mat(rows[row], cols[col]);
}

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteDenseLinearOperator<ValueType>::domain() const
{
    return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteDenseLinearOperator<ValueType>::range() const
{
    return m_rangeSpace;
}

template <typename ValueType>
bool DiscreteDenseLinearOperator<ValueType>::opSupportedImpl(
        Thyra::EOpTransp M_trans) const
{
    return (M_trans == Thyra::NOTRANS || M_trans == Thyra::TRANS
            || M_trans == Thyra::CONJ || M_trans == Thyra::CONJTRANS);
}

template <typename ValueType>
void DiscreteDenseLinearOperator<ValueType>::applyImpl(
        const Thyra::EOpTransp M_trans,
        const Thyra::MultiVectorBase<ValueType> &X_in,
        const Teuchos::Ptr<Thyra::MultiVectorBase<ValueType> > &Y_inout,
        const ValueType alpha,
        const ValueType beta) const
{
    typedef Thyra::Ordinal Ordinal;

    // Note: the name is VERY misleading: these asserts don't disappear in
    // release runs, and in case of failure throw exceptions rather than
    // abort.
    TEUCHOS_ASSERT(this->opSupported(M_trans));
    TEUCHOS_ASSERT(X_in.range()->isCompatible(*this->domain()));
    TEUCHOS_ASSERT(Y_inout->range()->isCompatible(*this->range()));
    TEUCHOS_ASSERT(Y_inout->domain()->isCompatible(*X_in.domain()));

    const Ordinal colCount = X_in.domain()->dim();

    // Loop over the input columns

    for (Ordinal col = 0; col < colCount; ++col) {
        // Get access the the elements of X_in's and Y_inout's column #col
        Thyra::ConstDetachedSpmdVectorView<ValueType> xVec(X_in.col(col));
        Thyra::DetachedSpmdVectorView<ValueType> yVec(Y_inout->col(col));
        const Teuchos::ArrayRCP<const ValueType> xArray(xVec.sv().values());
        const Teuchos::ArrayRCP<ValueType> yArray(yVec.sv().values());

        // Wrap the Trilinos array in an Armadillo vector. const_cast is used
        // because it's more natural to have a const arma::Col<ValueType> array
        // than an arma::Col<const ValueType> one.
        const arma::Col<ValueType> xCol(
                    const_cast<ValueType*>(xArray.get()), xArray.size(),
                    false /* copy_aux_mem */);
        arma::Col<ValueType> yCol(yArray.get(), yArray.size(), false);

        applyBuiltInImpl(static_cast<TranspositionMode>(M_trans),
                         xCol, yCol, alpha, beta);
    }
}
#endif // WITH_TRILINOS

template <typename ValueType>
void DiscreteDenseLinearOperator<ValueType>::applyBuiltInImpl(
        const TranspositionMode trans,
        const arma::Col<ValueType>& x_in,
        arma::Col<ValueType>& y_inout,
        const ValueType alpha,
        const ValueType beta) const
{
    if (columnCount() != x_in.n_rows && rowCount() != y_inout.n_rows)
        throw std::invalid_argument(
                "DiscreteDenseLinearOperator::applyBuiltInImpl(): "
                "incorrect vector length");

    if (beta == static_cast<ValueType>(0.))
        y_inout.fill(static_cast<ValueType>(0.));
    else
        y_inout *= beta;

    switch (trans)
    {
    case NO_TRANSPOSE:
        y_inout += alpha * m_mat * x_in;
        break;
    case CONJUGATE:
        y_inout += alpha * arma::conj(m_mat) * x_in;
        break;
    case TRANSPOSE:
        y_inout += alpha * m_mat.st() * x_in;
        break;
    case CONJUGATE_TRANSPOSE:
        y_inout += alpha * m_mat.t() * x_in;
        break;
    default:
        throw std::invalid_argument(
                "AcaApproximateLuInverse::applyBuiltInImpl(): "
                "invalid transposition mode");
    }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteDenseLinearOperator);

} // namespace Bempp
