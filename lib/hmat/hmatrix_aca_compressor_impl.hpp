// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_ACA_COMPRESSOR_IMPL_HPP
#define HMAT_HMATRIX_ACA_COMPRESSOR_IMPL_HPP

#include "hmatrix_aca_compressor.hpp"
#include "hmatrix_low_rank_data.hpp"
#include "scalar_traits.hpp"
#include <random>
#include <complex>
#include <cmath>

namespace hmat {

template <typename ValueType, int N>
void HMatrixAcaCompressor<ValueType, N>::compressBlock(
    const BlockClusterTreeNode<N> &blockClusterTreeNode,
    shared_ptr<HMatrixData<ValueType>> &hMatrixData) const {

  if (!blockClusterTreeNode.data().admissible) {
    m_hMatrixDenseCompressor.compressBlock(blockClusterTreeNode, hMatrixData);
    return;
  }

  std::vector<shared_ptr<arma::Mat<ValueType>>> previousColumns;
  std::vector<shared_ptr<arma::Mat<ValueType>>> previousRows;

  IndexRangeType rowClusterRange;
  IndexRangeType columnClusterRange;
  std::size_t numberOfRows;
  std::size_t numberOfColumns;

  getBlockClusterTreeNodeDimensions(blockClusterTreeNode, rowClusterRange,
                                    columnClusterRange, numberOfRows,
                                    numberOfColumns);

  IndexRangeType rowIndexRange;
  IndexRangeType columnIndexRange;

  for (int i = 0; i < m_maxRank; ++i) {

    std::size_t row = intRand(rowClusterRange);
    shared_ptr<arma::Mat<ValueType>> newCol(new arma::Mat<ValueType>());
    shared_ptr<arma::Mat<ValueType>> newRow(new arma::Mat<ValueType>());

    // Compute complete row
    rowIndexRange = {{row, row + 1}};
    columnIndexRange = columnClusterRange;

    evaluateMatMinusLowRank(blockClusterTreeNode, rowIndexRange,
                            columnIndexRange, *newRow, previousColumns,
                            previousRows);

    arma::uword maxRowInd;
    arma::uword maxColInd;

    arma::Mat<typename ScalarTraits<ValueType>::RealType> absMat = arma::abs(*newRow);

    if (absMat.max(maxRowInd,maxColInd) < 1E-12)
      continue;

    auto pivot = (*newRow)(row, maxColInd);

    *newRow = *newRow / (pivot);

    // Now evaluate column

    rowIndexRange = rowClusterRange;
    columnClusterRange = {{maxColInd, maxColInd + 1}};

    evaluateMatMinusLowRank(blockClusterTreeNode, rowIndexRange,
                            columnIndexRange, *newCol, previousColumns,
                            previousRows);

    previousRows.push_back(newRow);
    previousColumns.push_back(newCol);

    if (std::abs(pivot) < m_eps)
      break;
  }

  // Now fill up the hMatrixData structure

  hMatrixData.reset(new HMatrixLowRankData<ValueType>());

  arma::Mat<ValueType> &A =
      static_cast<HMatrixLowRankData<ValueType>*>(hMatrixData.get())->A();

  arma::Mat<ValueType> &B =
      static_cast<HMatrixLowRankData<ValueType>*>(hMatrixData.get())->B();

  A.resize(numberOfRows, previousColumns.size());
  B.resize(previousRows.size(), numberOfColumns);

  for (int i = 0; i < previousColumns.size(); ++i) {
    A.col(i) = *(previousColumns[i]);
    B.row(i) = *(previousRows[i]);
  }
}

template <typename ValueType, int N>
HMatrixAcaCompressor<ValueType, N>::HMatrixAcaCompressor(
    const DataAccessor<ValueType, N> &dataAccessor, double eps,
    unsigned int maxRank)
    : m_dataAccessor(dataAccessor), m_eps(eps), m_maxRank(maxRank),
      m_hMatrixDenseCompressor(dataAccessor) {}

template <typename ValueType, int N>
void HMatrixAcaCompressor<ValueType, N>::evaluateMatMinusLowRank(
    const BlockClusterTreeNode<N> &blockClusterTreeNode,
    const IndexRangeType &rowIndexRange, const IndexRangeType &columnIndexRange,
    arma::Mat<ValueType> &data,
    const std::vector<shared_ptr<arma::Mat<ValueType>>> &previousColumns,
    const std::vector<shared_ptr<arma::Mat<ValueType>>> &previousRows) const {

  auto rowClusterRange =
      blockClusterTreeNode.data().rowClusterTreeNode->data().indexRange;
  auto columnClusterRange =
      blockClusterTreeNode.data().columnClusterTreeNode->data().indexRange;

  m_dataAccessor.computeMatrixBlock(rowIndexRange, columnIndexRange,
                                    blockClusterTreeNode, data);

  auto rowStart = rowIndexRange[0] - rowClusterRange[0];
  auto rowEnd = rowIndexRange[1] - rowClusterRange[0];

  auto colStart = columnIndexRange[0] - columnClusterRange[0];
  auto colEnd = columnIndexRange[1] - columnIndexRange[0];

  for (int i = 0; i < previousColumns.size(); ++i)
    data -= (previousColumns[i]->rows(rowStart, rowEnd - 1)) *
            (previousRows[i]->cols(colStart, colEnd - 1));
}

template <typename ValueType, int N>
std::size_t
HMatrixAcaCompressor<ValueType, N>::intRand(const IndexRangeType &range) {
  static std::random_device generator;
  std::uniform_int_distribution<std::size_t> distribution(range[0],
                                                          range[1] - 1);
  return distribution(generator);
}
}

#endif
