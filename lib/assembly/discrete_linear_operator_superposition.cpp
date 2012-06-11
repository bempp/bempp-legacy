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

#include "discrete_linear_operator_superposition.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../common/not_implemented_error.hpp"

#ifdef WITH_TRILINOS
#include <Thyra_DetachedSpmdVectorView.hpp>
#include <Thyra_SpmdVectorSpaceDefaultBase.hpp>
#endif

namespace Bempp
{

template <typename ValueType>
DiscreteLinearOperatorSuperposition<ValueType>::
DiscreteLinearOperatorSuperposition(
        boost::ptr_vector<TermType>& terms,
        const std::vector<ValueType>& weights)
{
    // Check that all terms have the same dimensions
    if (terms.size() > 1) {
        const size_t nRows = terms[0].rowCount();
        const size_t nCols = terms[0].columnCount();
        for (size_t i = 1; i < terms.size(); ++i)
            if (terms[i].rowCount() != nRows || terms[i].columnCount() != nCols)
                throw std::invalid_argument(
                        "DiscreteLinearOperatorSuperposition::"
                        "DiscreteLinearOperatorSuperposition(): "
                        "all terms must have the same dimensions");
    }

    m_terms.transfer(m_terms.end(), terms);
    m_weights.insert(m_weights.end(), weights.begin(), weights.end());
#ifdef WITH_TRILINOS
    if (!m_terms.empty()) {
        m_domainSpace = m_terms[0].domain();
        m_rangeSpace = m_terms[0].range();
    }
#endif
}

template <typename ValueType>
void DiscreteLinearOperatorSuperposition<ValueType>::dump() const
{
    throw NotImplementedError(
                "DiscreteLinearOperatorSuperposition::dump(): "
                "not implemented yet");
}

template <typename ValueType>
arma::Mat<ValueType>
DiscreteLinearOperatorSuperposition<ValueType>::asMatrix() const
{
    arma::Mat<ValueType> result;
    if (!m_terms.empty()) {
        result = m_terms[0].asMatrix();
        for (size_t i = 1; i < m_terms.size(); ++i)
            result += m_terms[i].asMatrix() * m_weights[i];
    }
    return result;
}

template <typename ValueType>
unsigned int
DiscreteLinearOperatorSuperposition<ValueType>::rowCount() const
{
    if (m_terms.empty())
        return 0;
    else
        return m_terms[0].rowCount();
}

template <typename ValueType>
unsigned int
DiscreteLinearOperatorSuperposition<ValueType>::columnCount() const
{
    if (m_terms.empty())
        return 0;
    else
        return m_terms[0].columnCount();
}

template <typename ValueType>
void DiscreteLinearOperatorSuperposition<ValueType>::addBlock(
        const std::vector<int>& rows,
        const std::vector<int>& cols,
        const ValueType alpha,
        arma::Mat<ValueType>& block) const
{
    for (size_t i = 0; i < m_terms.size(); ++i)
        m_terms[i].addBlock(rows, cols, alpha * m_weights[i], block);
}

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteLinearOperatorSuperposition<ValueType>::domain() const
{
    return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
DiscreteLinearOperatorSuperposition<ValueType>::range() const
{
    return m_rangeSpace;
}

template <typename ValueType>
bool DiscreteLinearOperatorSuperposition<ValueType>::opSupportedImpl(
        Thyra::EOpTransp M_trans) const
{
    // TODO: implement remaining variants (transpose & conjugate transpose)
    return (M_trans == Thyra::NOTRANS);
}

template <typename ValueType>
void DiscreteLinearOperatorSuperposition<ValueType>::applyImpl(
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
void DiscreteLinearOperatorSuperposition<ValueType>::
applyBuiltInImpl(const TranspositionMode trans,
                 const arma::Col<ValueType>& x_in,
                 arma::Col<ValueType>& y_inout,
                 const ValueType alpha,
                 const ValueType beta) const
{
    if (beta == static_cast<ValueType>(0.))
        y_inout.fill(static_cast<ValueType>(0.));
    else
        y_inout *= beta;

    for (size_t i = 0; i < m_terms.size(); ++i)
        m_terms[i].apply(trans, x_in, y_inout,
                         alpha * m_weights[i],
                         1. /* "+ beta * y_inout" has already been done */ );
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteLinearOperatorSuperposition);

} // namespace Bempp
