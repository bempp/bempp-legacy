// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_DENSE_DATA_IMPL_HPP
#define HMAT_HMATRIX_DENSE_DATA_IMPL_HPP

#include "common.hpp"

#include "hmatrix_dense_data.hpp"
#include "eigen_fwd.hpp"

namespace hmat {

template <typename ValueType>
void HMatrixDenseData<ValueType>::apply(const Matrix<ValueType> &X,
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
    Y = alpha * m_A * X + beta * Y;
  else if (trans == TransposeMode::TRANS)
    Y = alpha * m_A.transpose() * X + beta * Y;
  else if (trans == TransposeMode::CONJ)
    Y = alpha * m_A.conjugate() * X + beta * Y;
  else
    Y = alpha * m_A.adjoint() * X + beta * Y;
}

template <typename ValueType>
const Matrix<ValueType> &HMatrixDenseData<ValueType>::A() const {
  return m_A;
}

template <typename ValueType>
Matrix<ValueType> &HMatrixDenseData<ValueType>::A() {
  return m_A;
}

template <typename ValueType> int HMatrixDenseData<ValueType>::rows() const {
  return m_A.rows();
}

template <typename ValueType> int HMatrixDenseData<ValueType>::cols() const {
  return m_A.cols();
}

template <typename ValueType> int HMatrixDenseData<ValueType>::rank() const {
  return m_A.cols();
}

template <typename ValueType>
DataBlockType HMatrixDenseData<ValueType>::type() const {

  return DENSE;

}

template <typename ValueType>
typename ScalarTraits<ValueType>::RealType
HMatrixDenseData<ValueType>::frobeniusNorm() const {
  return m_A.norm();
}

template <typename ValueType>
double HMatrixDenseData<ValueType>::memSizeKb() const {
  return sizeof(ValueType) * (this->rows()) * (this->cols()) / (1.0 * 1024);
}
}
#endif
