#ifndef HMAT_HMATRIX_COMPRESSOR_IMPL_HPP
#define HMAT_HMATRIX_COMPRESSOR_IMPL_HPP

#include "hmatrix_compressor.hpp"
#include "hmatrix_low_rank_data.hpp"

namespace hmat {

template <typename ValueType, int N>
HMatrixCompressor<ValueType, N>::HMatrixCompressor(double cutoff)
    : m_cutoff(cutoff) {}

template <typename ValueType, int N>
void HMatrixCompressor<ValueType, N>::compressBlock(
    const BlockClusterTreeNode<N> &blockClusterTreeNode,
    shared_ptr<HMatrixData<ValueType>> &hMatrixData) const {

  const auto &rowClusterData =
      blockClusterTreeNode.data().rowClusterTreeNode->data();
  const auto &colClusterData =
      blockClusterTreeNode.data().columnClusterTreeNode->data();

  double distance =
      rowClusterData.boundingBox.distance(colClusterData.boundingBox);

  if ((distance <= m_cutoff) || (!blockClusterTreeNode.data().admissible))
    compressBlockImpl(blockClusterTreeNode, hMatrixData);
  else {
    // If admissible create a zero HMatrix data block and return this
    auto numberOfRows =
        rowClusterData.indexRange[1] - rowClusterData.indexRange[0];
    auto numberOfCols =
        colClusterData.indexRange[1] - colClusterData.indexRange[0];
    hMatrixData.reset(new HMatrixLowRankData<ValueType>());
    static_cast<HMatrixLowRankData<ValueType> *>(hMatrixData.get())
        ->A()
        .resize(numberOfRows, 0);
    static_cast<HMatrixLowRankData<ValueType> *>(hMatrixData.get())
        ->B()
        .resize(numberOfCols, 0);
  }

  compressBlockImpl(blockClusterTreeNode, hMatrixData);
}
}
#endif
