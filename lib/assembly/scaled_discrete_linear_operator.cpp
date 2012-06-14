#include "scaled_discrete_linear_operator.hpp"

#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

template <typename ValueType>
ScaledDiscreteLinearOperator<ValueType>::ScaledDiscreteLinearOperator(
        ValueType multiplier, std::auto_ptr<Base>& op) :
    m_multiplier(multiplier), m_operator(op)
{
    if (!m_operator.get())
        throw std::invalid_argument(
            "ScaledDiscreteLinearOperator::ScaledDiscreteLinearOperator(): "
            "the wrapped operator must not be NULL");
}

template <typename ValueType>
arma::Mat<ValueType> ScaledDiscreteLinearOperator<ValueType>::asMatrix() const
{
    return m_multiplier * m_operator->asMatrix();
}

template <typename ValueType>
unsigned int ScaledDiscreteLinearOperator<ValueType>::rowCount() const
{
    return m_operator->rowCount();
}

template <typename ValueType>
unsigned int ScaledDiscreteLinearOperator<ValueType>::columnCount() const
{
    return m_operator->columnCount();
}

template <typename ValueType>
void ScaledDiscreteLinearOperator<ValueType>::addBlock(
        const std::vector<int>& rows,
        const std::vector<int>& cols,
        const ValueType alpha,
        arma::Mat<ValueType>& block) const
{
    m_operator->addBlock(rows, cols, m_multiplier * alpha, block);
}

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
ScaledDiscreteLinearOperator<ValueType>::domain() const
{
    return m_operator->domain();
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType> >
ScaledDiscreteLinearOperator<ValueType>::range() const
{
    return m_operator->range();
}

template <typename ValueType>
bool ScaledDiscreteLinearOperator<ValueType>::opSupportedImpl(
        Thyra::EOpTransp M_trans) const
{
    return m_operator->opSupported(M_trans);
}
#endif

template <typename ValueType>
void ScaledDiscreteLinearOperator<ValueType>::applyBuiltInImpl(
        const TranspositionMode trans,
        const arma::Col<ValueType>& x_in,
        arma::Col<ValueType>& y_inout,
        const ValueType alpha,
        const ValueType beta) const
{
    m_operator->apply(trans, x_in, y_inout,
                      m_multiplier * alpha, beta);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(ScaledDiscreteLinearOperator);

} // namespace Bempp
