// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_DENSE_DATA_HPP
#define HMAT_HMATRIX_DENSE_DATA_HPP

#include "common.hpp"
#include "hmatrix_data.hpp"
#include <armadillo>

namespace hmat {

template <typename ValueType> class HMatrixDenseData : HMatrixData<ValueType> {
public:
  void apply(const arma::Mat<ValueType> &X, arma::Mat<ValueType> &Y,
             TransposeMode trans, ValueType alpha, ValueType beta) const
      override;

  const arma::Mat<ValueType> &A() const;
  arma::Mat<ValueType> &A();

  int rows() const override;
  int cols() const override;
  int rank() const override;

  double memSizeKb() const override;

private:
  arma::Mat<ValueType> m_A;
};
}

#include "hmatrix_dense_data_impl.hpp"

#endif
