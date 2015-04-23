#include "scaled_discrete_boundary_operator.hpp"
#include "discrete_aca_boundary_operator.hpp"

#include "../common/complex_aux.hpp"
#include "../common/eigen_support.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp {

template <typename ValueType>
ScaledDiscreteBoundaryOperator<ValueType>::ScaledDiscreteBoundaryOperator(
    ValueType multiplier, const shared_ptr<const Base> &op)
    : m_multiplier(multiplier), m_operator(op) {
  if (!m_operator.get())
    throw std::invalid_argument(
        "ScaledDiscreteBoundaryOperator::ScaledDiscreteBoundaryOperator(): "
        "the wrapped operator must not be NULL");
}

template <typename ValueType>
Matrix<ValueType>
ScaledDiscreteBoundaryOperator<ValueType>::asMatrix() const {
  return m_multiplier * m_operator->asMatrix();
}

template <typename ValueType>
unsigned int ScaledDiscreteBoundaryOperator<ValueType>::rowCount() const {
  return m_operator->rowCount();
}

template <typename ValueType>
unsigned int ScaledDiscreteBoundaryOperator<ValueType>::columnCount() const {
  return m_operator->columnCount();
}

template <typename ValueType>
void ScaledDiscreteBoundaryOperator<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, Matrix<ValueType> &block) const {
  m_operator->addBlock(rows, cols, m_multiplier * alpha, block);
}

#ifdef WITH_AHMED
template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>>
ScaledDiscreteBoundaryOperator<ValueType>::asDiscreteAcaBoundaryOperator(
    double eps, int maximumRank, bool interleave) const {
  shared_ptr<const DiscreteBoundaryOperator<ValueType>> acaOp =
      m_operator->asDiscreteAcaBoundaryOperator(eps, maximumRank, interleave);
  return scaledAcaOperator(acaOp, m_multiplier);
}
#endif // WITH_AHMED

template <typename ValueType>
void ScaledDiscreteBoundaryOperator<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const Eigen::Ref<Vector<ValueType>> &x_in,
    Eigen::Ref<Vector<ValueType>> y_inout, const ValueType alpha,
    const ValueType beta) const {
  ValueType multiplier = m_multiplier;
  if (trans == CONJUGATE || trans == CONJUGATE_TRANSPOSE)
    multiplier = conj(multiplier);
  m_operator->apply(trans, x_in, y_inout, multiplier * alpha, beta);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(ScaledDiscreteBoundaryOperator);

} // namespace Bempp
