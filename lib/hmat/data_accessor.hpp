// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_DATA_ACCESSOR_HPP
#define HMAT_DATA_ACCESSOR_HPP

#include "common.hpp"
#include <armadillo>

namespace hmat {

template <typename ValueType, int N> class DataAccessor {

  virtual void computeMatrixBlock(
      const IndexRangeType& rowIndexRange,
      const IndexRangeType& columnIndexRange,
      const BlockClusterTreeNode<N>& blockClusterTreeNode, 
      arma::Mat<ValueType> &data) const = 0;
};
}

#endif
