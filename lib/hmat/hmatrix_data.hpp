// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_DATA_HPP
#define HMAT_HMATRIX_DATA_HPP

#include "common.hpp"
#include "scalar_traits.hpp"
#include "eigen_fwd.hpp"

#include <vector>

namespace hmat {

template <typename ValueType> class HMatrixData {
public:
  virtual void apply(const Eigen::Ref<Matrix<ValueType>> &X,
                     Eigen::Ref<Matrix<ValueType>> Y, TransposeMode trans,
                     ValueType alpha, ValueType beta) const = 0;

  virtual int rows() const = 0;
  virtual int cols() const = 0;
  virtual int rank() const = 0;

  virtual typename ScalarTraits<ValueType>::RealType frobeniusNorm() const = 0;

  virtual double memSizeKb() const = 0;

  virtual int numberOfElements() const = 0;

  virtual DataBlockType type() const = 0;
};
}

#endif
