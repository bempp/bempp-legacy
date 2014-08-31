// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_COMPRESSOR_HPP
#define HMAT_HMATRIX_COMPRESSOR_HPP

#include "common.hpp"
#include "block_cluster_tree.hpp"
#include <armadillo>

namespace hmat {

template <typename ValueType> class HMatrixData;

template <typename ValueType, int N> class HMatrixCompressor {
public:
  virtual void compressBlock(
      const BlockClusterTreeNode<N> &blockClusterTreeNode,
      shared_ptr<HMatrixData<ValueType>> &hMatrixData) const = 0;
};
}

#endif
