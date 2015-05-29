// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_DATA_IMPL_HPP
#define HMAT_HMATRIX_DATA_IMPL_HPP

#include "hmatrix_low_rank_data.hpp"
#include "eigen_fwd.hpp"

namespace hmat {

template <typename ValueType>
const Matrix<ValueType> &HMatrixLowRankData<ValueType>::A() const {
  return m_A;
}

template <typename ValueType>
Matrix<ValueType> &HMatrixLowRankData<ValueType>::A() {
  return m_A;
}

template <typename ValueType>
const Matrix<ValueType> &HMatrixLowRankData<ValueType>::B() const {
  return m_B;
}

template <typename ValueType>
Matrix<ValueType> &HMatrixLowRankData<ValueType>::B() {
  return m_B;
}

template <typename ValueType> int HMatrixLowRankData<ValueType>::rows() const {
  return m_A.rows();
}

template <typename ValueType> int HMatrixLowRankData<ValueType>::cols() const {

  return m_B.cols();
}

template <typename ValueType> int HMatrixLowRankData<ValueType>::rank() const {

  return m_A.cols();
}

template <typename ValueType>
typename ScalarTraits<ValueType>::RealType
HMatrixLowRankData<ValueType>::frobeniusNorm() const {

  auto aHa = m_A.adjoint() * m_A;

  Matrix<ValueType> result(1, 1);

  for (int i = 0; i < m_B.cols(); ++i) {
    auto col = m_B.col(i);
    result += col.adjoint() * aHa * col;
  }

  return std::sqrt(std::real(result(0, 0)));
}

template <typename ValueType>
double HMatrixLowRankData<ValueType>::memSizeKb() const {

  return sizeof(ValueType) * (this->rows() + this->cols()) * this->rank() /
         (1.0 * 1024);
}

template <typename ValueType>
void HMatrixLowRankData<ValueType>::apply(const Matrix<ValueType> &X,
                                          Matrix<ValueType> &Y,
                                          TransposeMode trans, ValueType alpha,
                                          ValueType beta) const {

  if (beta == ValueType(0))
    Y.setZero();
  if (alpha == ValueType(0)) {
    Y = beta * Y;
    return;
  }

  if (trans == TransposeMode::NOTRANS)
    Y = alpha * m_A * (m_B * X) + beta * Y;
  else if (trans == TransposeMode::TRANS)
    Y = alpha * m_B.transpose() * (m_A.transpose() * X) + beta * Y;
  else if (trans == TransposeMode::CONJ)
    Y = alpha * m_A.conjugate() * (m_B.conjugate()* X) + beta * Y;
  else
    Y = alpha * m_B.adjoint() * (m_A.adjoint() * X) + beta * Y;
}
}

#endif
