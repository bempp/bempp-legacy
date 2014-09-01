// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_DATA_IMPL_HPP
#define HMAT_HMATRIX_DATA_IMPL_HPP

#include "hmatrix_low_rank_data.hpp"

namespace hmat {

template <typename ValueType>
const arma::Mat<ValueType> &HMatrixLowRankData<ValueType>::A() const {
  return m_A;
}

template <typename ValueType>
arma::Mat<ValueType> &HMatrixLowRankData<ValueType>::A() {
  return m_A;
}

template <typename ValueType>
const arma::Mat<ValueType> &HMatrixLowRankData<ValueType>::B() const {
  return m_B;
}

template <typename ValueType>
arma::Mat<ValueType> &HMatrixLowRankData<ValueType>::B() {
  return m_B;
}

template <typename ValueType> int HMatrixLowRankData<ValueType>::rows() const {
  return m_A.n_rows;
}

template <typename ValueType> int HMatrixLowRankData<ValueType>::cols() const {

  return m_B.n_cols;
}

template <typename ValueType> int HMatrixLowRankData<ValueType>::rank() const {

  return m_A.n_cols;
}

template <typename ValueType>
double HMatrixLowRankData<ValueType>::memSizeKb() const {

  return sizeof(ValueType) * (this->rows() + this->cols()) * this->rank() /
         (1.0 * 1024);
}

template <typename ValueType>
void HMatrixLowRankData<ValueType>::apply(const arma::Mat<ValueType> &X,
                                          arma::Mat<ValueType> &Y,
                                          TransposeMode trans, ValueType alpha,
                                          ValueType beta) const {

  arma::subview<ValueType> xsub = X.submat(arma::span::all, arma::span::all);
  arma::subview<ValueType> ysub = Y.submat(arma::span::all, arma::span::all);

  this->apply(xsub, ysub, trans, alpha, beta);
}
template <typename ValueType>
void HMatrixLowRankData<ValueType>::apply(const arma::subview<ValueType> &X,
                                          arma::subview<ValueType> &Y,
                                          TransposeMode trans, ValueType alpha,
                                          ValueType beta) const {
  if (beta == ValueType(0))
    Y.zeros();
  if (alpha == ValueType(0)) {
    Y = beta * Y;
    return;
  }

  if (trans == TransposeMode::NOTRANS)
    Y = alpha * m_A * (m_B * X) + beta * Y;
  else if (trans == TransposeMode::TRANS)
    Y = alpha * m_B.st() * (m_A.st() * X) + beta * Y;
  else if (trans == TransposeMode::CONJ)
    Y = alpha * arma::conj(m_A) * (arma::conj(m_B) * X) + beta * Y;
  else
    Y = alpha * m_B.t() * (m_A.t() * X) + beta * Y;
}
}

#endif
