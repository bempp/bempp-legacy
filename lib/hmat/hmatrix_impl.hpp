// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_IMPL_HPP
#define HMAT_HMATRIX_IMPL_HPP

#include "hmatrix.hpp"
#include "hmatrix_data.hpp"

namespace hmat {

template <typename ValueType, int N>
HMatrix<ValueType, N>::HMatrix(
    const shared_ptr<BlockClusterTree<N>> &blockClusterTree)
    : m_blockClusterTree(blockClusterTree) {}

template <typename ValueType, int N>
HMatrix<ValueType, N>::HMatrix(
    const shared_ptr<BlockClusterTree<N>> &blockClusterTree,
    const HMatrixCompressor<ValueType, N> &hMatrixCompressor,
    const DataAccessor<ValueType, N> &dataAccessor)
    : HMatrix<ValueType, N>(blockClusterTree) {
  initialize(hMatrixCompressor, dataAccessor);
}

template <typename ValueType, int N>
void HMatrix<ValueType, N>::initialize(
    const HMatrixCompressor<ValueType, N> &hMatrixCompressor,
    const DataAccessor<ValueType, N> &dataAccessor) {

  reset();

  for (auto node: m_blockClusterTree->leafNodes()){
    std::shared_ptr<HMatrixData<ValueType>> nodeData;
    hMatrixCompressor.compressBlock(*node,nodeData);
    m_hMatrixData[node] = nodeData;
  }

}
template <typename ValueType, int N> void HMatrix<ValueType, N>::reset() {
  m_hMatrixData.clear();
}
}

#endif
