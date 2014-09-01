// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_COMPRESSED_MATRIX_HPP
#define HMAT_COMPRESSED_MATRIX_HPP

#include "common.hpp"
#include <armadillo>

namespace hmat {

template <typename ValueType> class CompressedMatrix {
public:
  virtual std::size_t rows() const = 0;
  virtual std::size_t columns() const = 0;

  virtual arma::Mat<ValueType>
  permuteMatToHMatDofs(const arma::Mat<ValueType> &mat,
                   RowColSelector rowOrColumn) const = 0;

  virtual arma::Mat<ValueType>
  permuteMatToOriginalDofs(const arma::Mat<ValueType> &mat,
                       RowColSelector rowOrColumn) const = 0;

  virtual void apply(const arma::Mat<ValueType> &X, arma::Mat<ValueType> &Y,
             TransposeMode trans, ValueType alpha, ValueType beta) const;
};
}

#endif

