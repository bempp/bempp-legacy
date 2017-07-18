// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_COMPRESSOR_HPP
#define HMAT_HMATRIX_COMPRESSOR_HPP

#include "block_cluster_tree.hpp"
#include "common.hpp"
#include "eigen_fwd.hpp"

namespace hmat {

template <typename ValueType> class HMatrixData;

template <typename ValueType, int N> class HMatrixCompressor {
public:
  HMatrixCompressor(double cutoff);
  void compressBlock(const BlockClusterTreeNode<N> &blockClusterTreeNode,
                     shared_ptr<HMatrixData<ValueType>> &hMatrixData) const;

protected:
  virtual void
  compressBlockImpl(const BlockClusterTreeNode<N> &blockClusterTreeNode,
                    shared_ptr<HMatrixData<ValueType>> &hMatrixData) const = 0;

private:
  double m_cutoff;
};
}

#include "hmatrix_compressor_impl.hpp"

#endif
