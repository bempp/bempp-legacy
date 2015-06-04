// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_DENSE_DATA_HPP
#define HMAT_HMATRIX_DENSE_DATA_HPP

#include "common.hpp"
#include "hmatrix_data.hpp"
#include "eigen_fwd.hpp"

namespace hmat {

template <typename ValueType>
class HMatrixDenseData : public HMatrixData<ValueType> {
public:
  void apply(const Eigen::Ref<Matrix<ValueType>> &X,
             Eigen::Ref<Matrix<ValueType>> Y, TransposeMode trans,
             ValueType alpha, ValueType beta) const override;

  const Matrix<ValueType> &A() const;
  Matrix<ValueType> &A();

  int rows() const override;
  int cols() const override;
  int rank() const override;

  typename ScalarTraits<ValueType>::RealType frobeniusNorm() const override;

  double memSizeKb() const override;

  DataBlockType type() const override;

private:
  Matrix<ValueType> m_A;
};
}

#include "hmatrix_dense_data_impl.hpp"

#endif
