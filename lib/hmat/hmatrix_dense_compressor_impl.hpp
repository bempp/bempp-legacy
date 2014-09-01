// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_DENSE_COMPRESSOR_IMPL_HPP
#define HMAT_HMATRIX_DENSE_COMPRESSOR_IMPL_HPP

#include "hmatrix_dense_compressor.hpp"
#include "hmatrix_dense_data.hpp"

namespace hmat {

template <typename ValueType, int N>
HMatrixDenseCompressor<ValueType, N>::HMatrixDenseCompressor(
    const DataAccessor<ValueType, N> &dataAccessor)
    : m_dataAccessor(dataAccessor) {}

template <typename ValueType, int N>
void HMatrixDenseCompressor<ValueType, N>::compressBlock(
    const BlockClusterTreeNode<N> &blockClusterTreeNode,
    shared_ptr<HMatrixData<ValueType>> &hMatrixData) const {

  auto rowIndexRange =
      blockClusterTreeNode.data().rowClusterTreeNode->data().indexRange;

  auto columnIndexRange =
      blockClusterTreeNode.data().columnClusterTreeNode->data().indexRange;

  hMatrixData.reset(new HMatrixDenseData<ValueType>());

  arma::Mat<ValueType> &A =
      static_cast<HMatrixDenseData<ValueType> *>(hMatrixData.get())->A();

  m_dataAccessor.computeMatrixBlock(rowIndexRange, columnIndexRange,
                                    blockClusterTreeNode, A);
}
}

#endif
