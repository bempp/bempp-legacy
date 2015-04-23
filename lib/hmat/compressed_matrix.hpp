// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_COMPRESSED_MATRIX_HPP
#define HMAT_COMPRESSED_MATRIX_HPP

#include "common.hpp"
#include "eigen_fwd.hpp"

namespace hmat {

template <typename ValueType> class CompressedMatrix {
public:
  virtual std::size_t rows() const = 0;
  virtual std::size_t columns() const = 0;

  virtual Matrix<ValueType>
  permuteMatToHMatDofs(const Matrix<ValueType> &mat,
                       RowColSelector rowOrColumn) const = 0;

  virtual Matrix<ValueType>
  permuteMatToOriginalDofs(const Matrix<ValueType> &mat,
                           RowColSelector rowOrColumn) const = 0;

  virtual void apply(const Matrix<ValueType> &X, Matrix<ValueType> &Y,
                     TransposeMode trans, ValueType alpha,
                     ValueType beta) const = 0;
};
}

#endif
