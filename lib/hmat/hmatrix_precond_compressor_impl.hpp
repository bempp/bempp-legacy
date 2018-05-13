// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_PRECOND_COMPRESSOR_IMPL_HPP
#define HMAT_HMATRIX_PRECOND_COMPRESSOR_IMPL_HPP

#include "hmatrix_precond_compressor.hpp"

namespace hmat {

template <typename ValueType, int N>
HMatrixPrecondCompressor<ValueType, N>::HMatrixPrecondCompressor(
    const DataAccessor<ValueType, N> &dataAccessor)
    : m_dataAccessor(dataAccessor) {}

template <typename ValueType, int N>
void HMatrixPrecondCompressor<ValueType, N>::compressBlock(
    const BlockClusterTreeNode<N> &blockClusterTreeNode,
    shared_ptr<HMatrixData<ValueType>> &hMatrixData) const {

  auto rowIndexRange =
      blockClusterTreeNode.data().rowClusterTreeNode->data().indexRange;

  auto columnIndexRange =
      blockClusterTreeNode.data().columnClusterTreeNode->data().indexRange;

  auto numberOfRowIndices = rowIndexRange[1] - rowIndexRange[0];
  auto numberOfColIndices = colIndexRange[1] - colIndexRange[0];

  hMatrixData.reset(new HMatrixLowRankData<ValueType>());

  Matrix<ValueType> &A =
      static_cast<HMatrixLowRankData<ValueType> *>(hMatrixData.get())->A();

  Matrix<ValueType> &B =
      static_cast<HMatrixLowRankData<ValueType> *>(hMatrixData.get())->B();

  Matrix<ValueType> row, col;

  const size_t maxRowIndices = 5;
  const size_t maxColumnIndices = 5;

  size_t rowOffset = maxRowIndices;
  size_t columnOffset = maxColumnIndices;

  if (maxRowIndices > numberOfRowIndices)
    rowOffset = numberOfRowIndices;
  if (maxColumnIndices > numberOfColumnIndices)
    columnOffset = numberOfColumnIndices;

  // Select the first few rows or columns
  IndexRangeType selectedRowIndexRange =
      {{rowIndexRange[0], rowIndexRange[0] + rowOffset}};
  IndexRangeType selectedColumnIndexRange =
      {{columnIndexRange[0], columnIndexRange[0] + columnOffset}};

  m_dataAccessor.computeMatrixBlock(selectedRowIndexRange, columnIndexRange,
                                    blockClusterTreeNode, row);

  m_dataAccessor.computeMatrixBlock(rowIndexRange, selectedColumnIndexRange,
                                    blockClusterTreeNode, col);
}
}

#endif
