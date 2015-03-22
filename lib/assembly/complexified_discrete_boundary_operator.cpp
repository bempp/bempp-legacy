#include "complexified_discrete_boundary_operator.hpp"
#include "discrete_aca_boundary_operator.hpp"
#include "../common/eigen_support.hpp"

#include "../fiber/explicit_instantiation.hpp"

#ifdef WITH_TRILINOS
#include <Thyra_DefaultSpmdVectorSpace_decl.hpp>
#endif

namespace Bempp {

template <typename RealType>
ComplexifiedDiscreteBoundaryOperator<RealType>::
    ComplexifiedDiscreteBoundaryOperator(
        const shared_ptr<const DiscreteBoundaryOperator<RealType>> &op)
    : m_operator(op)
#ifdef WITH_TRILINOS
      ,
      m_domainSpace(
          Thyra::defaultSpmdVectorSpace<ValueType>(op->columnCount())),
      m_rangeSpace(Thyra::defaultSpmdVectorSpace<ValueType>(op->rowCount()))
#endif
{
  if (!m_operator.get())
    throw std::invalid_argument("ComplexifiedDiscreteBoundaryOperator::"
                                "ComplexifiedDiscreteBoundaryOperator(): "
                                "the wrapped operator must not be NULL");
}

template <typename RealType>
Matrix<std::complex<RealType>>
ComplexifiedDiscreteBoundaryOperator<RealType>::asMatrix() const {
  Matrix<ValueType> result(rowCount(), columnCount());
  result.setZero();
  result.real() = m_operator->asMatrix();
  return result;
}

template <typename RealType>
unsigned int ComplexifiedDiscreteBoundaryOperator<RealType>::rowCount() const {
  return m_operator->rowCount();
}

template <typename RealType>
unsigned int
ComplexifiedDiscreteBoundaryOperator<RealType>::columnCount() const {
  return m_operator->columnCount();
}

template <typename RealType>
void ComplexifiedDiscreteBoundaryOperator<RealType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, Matrix<ValueType> &block) const {
  throw std::runtime_error("ComplexifiedDiscreteBoundaryOperator::addBlock(): "
                           "not implemented yet");
}

#ifdef WITH_TRILINOS
template <typename RealType>
Teuchos::RCP<const Thyra::VectorSpaceBase<std::complex<RealType>>>
ComplexifiedDiscreteBoundaryOperator<RealType>::domain() const {
  return m_domainSpace;
}

template <typename RealType>
Teuchos::RCP<const Thyra::VectorSpaceBase<std::complex<RealType>>>
ComplexifiedDiscreteBoundaryOperator<RealType>::range() const {
  return m_rangeSpace;
}

template <typename RealType>
bool ComplexifiedDiscreteBoundaryOperator<RealType>::opSupportedImpl(
    Thyra::EOpTransp M_trans) const {
  return m_operator->opSupported(M_trans);
}
#endif

template <typename RealType>
void ComplexifiedDiscreteBoundaryOperator<RealType>::applyBuiltInImpl(
    const TranspositionMode trans, const Vector<ValueType> &x_in,
    Vector<ValueType> &y_inout, const ValueType alpha,
    const ValueType beta) const {
  if (beta == static_cast<ValueType>(0.))
    y_inout.setZero();
  else
    y_inout *= beta;

  Vector<RealType> x_re = x_in.real();
  Vector<RealType> x_im = x_in.imag();

  Vector<RealType> y_re = y_inout.real();
  Vector<RealType> y_im = y_inout.imag();

  m_operator->apply(trans, x_re, y_re, alpha.real(), 1.);
  m_operator->apply(trans, x_im, y_re, -alpha.imag(), 1.);

  m_operator->apply(trans, x_re, y_im, alpha.imag(), 1.);
  m_operator->apply(trans, x_im, y_im, alpha.real(), 1.);

  y_inout.set_real(y_re);
  y_inout.set_imag(y_im);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_REAL_ONLY(
    ComplexifiedDiscreteBoundaryOperator);

} // namespace Bempp
