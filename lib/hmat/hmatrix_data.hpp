// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_DATA_HPP
#define HMAT_HMATRIX_DATA_HPP

namespace hmat {

template <typename ValueType> class HMatrixData {
public:
  virtual void apply(const arma::Mat<ValueType> &X, arma::Mat<ValueType> &Y,
                     TransposeMode trans, ValueType alpha,
                     ValueType beta) const = 0;

  virtual void apply(const arma::subview<ValueType> &X,
                     arma::subview<ValueType> &Y, TransposeMode trans,
                     ValueType alpha, ValueType beta) const = 0;

  virtual int rows() const = 0;
  virtual int cols() const = 0;
  virtual int rank() const = 0;

  virtual double memSizeKb() const = 0;
};
}

#endif
