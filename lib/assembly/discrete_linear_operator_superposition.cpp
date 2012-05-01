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

namespace Bempp
{

template <typename ValueType>
DiscreteLinearOperatorSuperposition<ValueType>::
DiscreteLinearOperatorSuperposition(
        boost::ptr_vector<TermType>& terms,
        const std::vector<ValueType>& multipliers)
{
    // Check that all terms have the same dimensions
    if (terms.size() > 1) {
        const int nRows = terms[0].rowCount();
        const int nCols = terms[0].columnCount();
        for (int i = 1; i < terms.size(); ++i)
            if (terms[i].rowCount() != nRows || terms[i].columnCount() != nCols)
                throw std::invalid_argument(
                        "DiscreteLinearOperatorSuperposition::"
                        "DiscreteLinearOperatorSuperposition(): "
                        "all terms must have the same dimensions");
    }

    m_terms.transfer(m_terms.end(), terms);
    m_multipliers.insert(m_multipliers.end(),multipliers.begin(),multipliers.end());
#ifdef WITH_TRILINOS
    if (!terms.empty()) {
        m_domainSpace = terms[0].domain();
        m_rangeSpace = terms[0].range();
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
        for (int i = 1; i < m_terms.size(); ++i)
            result += m_terms[i].asMatrix() * m_multipliers[i];
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
    for (int i = 0; i < m_terms.size(); ++i)
        m_terms[i].addBlock(rows, cols, alpha * m_multipliers[i], block);
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
    throw std::runtime_error("DiscreteLinearOperatorSuperposition::"
                             "applyImpl(): not implemented yet");
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

    for (int i = 0; i < m_terms.size(); ++i)
        m_terms[i].apply(trans, x_in, y_inout,
                         alpha * m_multipliers[i],
                         1. /* "+ beta * y_inout" has already been done */ );
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteLinearOperatorSuperposition);

} // namespace Bempp
