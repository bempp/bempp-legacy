#include "transposed_discrete_boundary_operator.hpp"
#include "discrete_aca_boundary_operator.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../common/eigen_support.hpp"

namespace Bempp {

template <typename ValueType>
TransposedDiscreteBoundaryOperator<ValueType>::
    TransposedDiscreteBoundaryOperator(TranspositionMode trans,
                                       const shared_ptr<const Base> &op)
    : m_trans(trans), m_operator(op) {
  if (!m_operator.get())
    throw std::invalid_argument("TransposedDiscreteBoundaryOperator::"
                                "TransposedDiscreteBoundaryOperator(): "
                                "the wrapped operator must not be NULL");
  if (m_trans != NO_TRANSPOSE && m_trans != TRANSPOSE && m_trans != CONJUGATE &&
      m_trans != CONJUGATE_TRANSPOSE)
    throw std::runtime_error("TransposedDiscreteBoundaryOperator::"
                             "TransposedDiscreteBoundaryOperator(): "
                             "invalid transposition mode");
}

template <typename ValueType>
Matrix<ValueType>
TransposedDiscreteBoundaryOperator<ValueType>::asMatrix() const {
  Matrix<ValueType> origMatrix = m_operator->asMatrix();
  switch (m_trans) {
  case NO_TRANSPOSE:
    return origMatrix;
  case TRANSPOSE:
    return origMatrix.transpose();
  case CONJUGATE:
    return origMatrix.conjugate();
  case CONJUGATE_TRANSPOSE:
    return origMatrix.adjoint();
  default:
    throw std::runtime_error("TransposedDiscreteBoundaryOperator::asMatrix(): "
                             "invalid transposition mode");
  }
}

template <typename ValueType>
bool TransposedDiscreteBoundaryOperator<ValueType>::isTransposed() const {
  return m_trans == TRANSPOSE || m_trans == CONJUGATE_TRANSPOSE;
}

template <typename ValueType>
unsigned int TransposedDiscreteBoundaryOperator<ValueType>::rowCount() const {
  return isTransposed() ? m_operator->columnCount() : m_operator->rowCount();
}

template <typename ValueType>
unsigned int
TransposedDiscreteBoundaryOperator<ValueType>::columnCount() const {
  return isTransposed() ? m_operator->rowCount() : m_operator->columnCount();
}

template <typename ValueType>
void TransposedDiscreteBoundaryOperator<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, Matrix<ValueType> &block) const {
  m_operator->addBlock(cols, rows, alpha, block);
}

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
TransposedDiscreteBoundaryOperator<ValueType>::domain() const {
  return isTransposed() ? m_operator->range() : m_operator->domain();
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
TransposedDiscreteBoundaryOperator<ValueType>::range() const {
  return isTransposed() ? m_operator->domain() : m_operator->range();
}

template <typename ValueType>
bool TransposedDiscreteBoundaryOperator<ValueType>::opSupportedImpl(
    Thyra::EOpTransp M_trans) const {
  // Bitwise xor. We use the fact that bit 0 of M_trans denotes
  // conjugation, and bit 1 -- transposition.
  return m_operator->opSupported(Thyra::EOpTransp(M_trans ^ m_trans));
}
#endif

template <typename ValueType>
void TransposedDiscreteBoundaryOperator<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const Vector<ValueType> &x_in,
    Vector<ValueType> &y_inout, const ValueType alpha,
    const ValueType beta) const {
  // Bitwise xor. We use the fact that bit 0 of M_trans denotes
  // conjugation, and bit 1 -- transposition.
  m_operator->apply(TranspositionMode(trans ^ m_trans), x_in, y_inout, alpha,
                    beta);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(TransposedDiscreteBoundaryOperator);

} // namespace Bempp
