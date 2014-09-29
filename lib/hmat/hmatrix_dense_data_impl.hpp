// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_DENSE_DATA_IMPL_HPP
#define HMAT_HMATRIX_DENSE_DATA_IMPL_HPP

#include "common.hpp"
#include <armadillo>

#include "hmatrix_dense_data.hpp"

namespace hmat {

template <typename ValueType>
void HMatrixDenseData<ValueType>::apply(const arma::Mat<ValueType> &X,
                                        arma::Mat<ValueType> &Y,
                                        TransposeMode trans, ValueType alpha,
                                        ValueType beta) const {

  arma::subview<ValueType> xsub = X.submat(arma::span::all, arma::span::all);
  arma::subview<ValueType> ysub = Y.submat(arma::span::all, arma::span::all);

  this->apply(xsub, ysub, trans, alpha, beta);
}

template <typename ValueType>
void HMatrixDenseData<ValueType>::apply(const arma::subview<ValueType> &X,
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
    Y = alpha * m_A * X + beta * Y;
  else if (trans == TransposeMode::TRANS)
    Y = alpha * m_A.st() * X + beta * Y;
  else if (trans == TransposeMode::CONJ)
    Y = alpha * arma::conj(m_A) * X + beta * Y;
  else
    Y = alpha * m_A.t() * X + beta * Y;
}

template <typename ValueType>
const arma::Mat<ValueType> &HMatrixDenseData<ValueType>::A() const {
  return m_A;
}

template <typename ValueType>
arma::Mat<ValueType> &HMatrixDenseData<ValueType>::A() {
  return m_A;
}

template <typename ValueType> int HMatrixDenseData<ValueType>::rows() const {
  return m_A.n_rows;
}

template <typename ValueType> int HMatrixDenseData<ValueType>::cols() const {
  return m_A.n_cols;
}

template <typename ValueType> int HMatrixDenseData<ValueType>::rank() const {
  return m_A.n_cols;
}

template <typename ValueType>
typename ScalarTraits<ValueType>::RealType HMatrixDenseData<ValueType>::frobeniusNorm() const {
  return norm(m_A,"fro");
}


template <typename ValueType>
double HMatrixDenseData<ValueType>::memSizeKb() const {
  return sizeof(ValueType) * (this->rows()) * (this->cols()) / (1.0 * 1024);
}
}
#endif
