// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_LOW_RANK_DATA_HPP
#define HMAT_HMATRIX_LOW_RANK_DATA_HPP

#include "common.hpp"
#include "hmatrix_data.hpp"
#include <armadillo>

namespace hmat {

template <typename ValueType>
class HMatrixLowRankData : public HMatrixData<ValueType> {

public:
  void apply(const arma::Mat<ValueType> &X, arma::Mat<ValueType> &Y,
             TransposeMode trans, ValueType alpha, ValueType beta) const
      override;

  void apply(const arma::subview<ValueType> &X, arma::subview<ValueType> &Y,
             TransposeMode trans, ValueType alpha, ValueType beta) const
      override;

  const arma::Mat<ValueType> &A() const;
  arma::Mat<ValueType> &A();

  const arma::Mat<ValueType> &B() const;
  arma::Mat<ValueType> &B();

  int rows() const override;
  int cols() const override;
  int rank() const override;

  double memSizeKb() const override;

private:
  arma::Mat<ValueType> m_A;
  arma::Mat<ValueType> m_B;
};
}

#include "hmatrix_low_rank_data_impl.hpp"

#endif
