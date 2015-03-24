// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_ACA_COMPRESSOR_IMPL_HPP
#define HMAT_HMATRIX_ACA_COMPRESSOR_IMPL_HPP

#include "hmatrix_aca_compressor.hpp"
#include "hmatrix_low_rank_data.hpp"
#include "eigen_fwd.hpp"
#include "scalar_traits.hpp"
#include <random>
#include <complex>
#include <cmath>
#include <algorithm>

namespace hmat {

template <typename ValueType, int N>
void HMatrixAcaCompressor<ValueType, N>::compressBlock(
    const BlockClusterTreeNode<N> &blockClusterTreeNode,
    shared_ptr<HMatrixData<ValueType>> &hMatrixData) const {

  if (!blockClusterTreeNode.data().admissible) {
    m_hMatrixDenseCompressor.compressBlock(blockClusterTreeNode, hMatrixData);
    return;
  }

  IndexRangeType rowClusterRange;
  IndexRangeType columnClusterRange;
  std::size_t numberOfRows;
  std::size_t numberOfColumns;

  getBlockClusterTreeNodeDimensions(blockClusterTreeNode, rowClusterRange,
                                    columnClusterRange, numberOfRows,
                                    numberOfColumns);

  hMatrixData.reset(new HMatrixLowRankData<ValueType>());

  Matrix<ValueType> &A =
      static_cast<HMatrixLowRankData<ValueType> *>(hMatrixData.get())->A();

  Matrix<ValueType> &B =
      static_cast<HMatrixLowRankData<ValueType> *>(hMatrixData.get())->B();

  A.resize(numberOfRows, m_resizeThreshold);
  B.resize(m_resizeThreshold, numberOfColumns);

  std::set<std::size_t> previousRowIndices;
  std::set<std::size_t> previousColumnIndices;

  std::size_t iterationLimit =
      std::min(static_cast<std::size_t>(m_maxRank),
               std::min(numberOfRows, numberOfColumns));

  std::size_t rankCount = 0;

  std::size_t sizeMultiplier = 0;

  for (int i = 0; i < iterationLimit; ++i) {

    std::size_t row = randomIndex(rowClusterRange, previousRowIndices);

    // Compute complete row

    Matrix<ValueType> newRow;
    Matrix<ValueType> newCol;

    IndexRangeType rowIndexRange = {{row, row + 1}};
    IndexRangeType columnIndexRange = columnClusterRange;

    evaluateMatMinusLowRank(blockClusterTreeNode, rowIndexRange,
                            columnIndexRange, newRow, A, B);

    std::ptrdiff_t maxRowInd;
    std::ptrdiff_t maxColInd;

    Matrix<typename ScalarTraits<ValueType>::RealType> absMat =
        newRow.cwiseAbs();

    if (absMat.maxCoeff(maxRowInd, maxColInd) < 1E-12)
      continue; // Row is effectively zero

    auto pivot = newRow(0, maxColInd);

    newRow = newRow / (pivot);

    // Now evaluate column

    maxColInd += columnClusterRange[0]; // Map back to original variables
    rowIndexRange = rowClusterRange;
    columnIndexRange = {{maxColInd, maxColInd + 1}};

    evaluateMatMinusLowRank(blockClusterTreeNode, rowIndexRange,
                            columnIndexRange, newCol, A, B);

    auto frobeniousNorm = hMatrixData->frobeniusNorm();

    if (rankCount == A.cols()) {
      sizeMultiplier++;
      A.insert_cols(sizeMultiplier * m_resizeThreshold, m_resizeThreshold);
      B.insert_rows(sizeMultiplier * m_resizeThreshold, m_resizeThreshold);
    }

    A.col(rankCount) = newCol;
    B.row(rankCount) = newRow;

    rankCount++;

    if (newCol.norm() * newRow.norm() < m_eps * frobeniousNorm)
      break;
  }
  if (A.cols() - rankCount > 0) {
      A = A.block(0,0,A.rows(),rankCount);
      B = B.block(0,0,B.rows(),rankCount);
  }
}

template <typename ValueType, int N>
HMatrixAcaCompressor<ValueType, N>::HMatrixAcaCompressor(
    const DataAccessor<ValueType, N> &dataAccessor, double eps,
    unsigned int maxRank, unsigned int resizeThreshold)
    : m_dataAccessor(dataAccessor), m_eps(eps), m_maxRank(maxRank),
      m_resizeThreshold(resizeThreshold),
      m_hMatrixDenseCompressor(dataAccessor) {}

template <typename ValueType, int N>
void HMatrixAcaCompressor<ValueType, N>::evaluateMatMinusLowRank(
    const BlockClusterTreeNode<N> &blockClusterTreeNode,
    const IndexRangeType &rowIndexRange, const IndexRangeType &columnIndexRange,
    Matrix<ValueType> &data, const Matrix<ValueType> &A,
    const Matrix<ValueType> &B) const {

  auto rowClusterRange =
      blockClusterTreeNode.data().rowClusterTreeNode->data().indexRange;
  auto columnClusterRange =
      blockClusterTreeNode.data().columnClusterTreeNode->data().indexRange;

  m_dataAccessor.computeMatrixBlock(rowIndexRange, columnIndexRange,
                                    blockClusterTreeNode, data);

  auto rowStart = rowIndexRange[0] - rowClusterRange[0];
  auto rowEnd = rowIndexRange[1] - rowClusterRange[0];

  auto colStart = columnIndexRange[0] - columnClusterRange[0];
  auto colEnd = columnIndexRange[1] - columnClusterRange[0];

  data = data - A.block(rowStart, 0, rowEnd - rowStart, A.cols()) *
                    B.submat(0, colStart, B.rows(), colEnd - colStart);
}

template <typename ValueType, int N>
std::size_t HMatrixAcaCompressor<ValueType, N>::randomIndex(
    const IndexRangeType &range, std::set<std::size_t> &previousIndices) {

  std::size_t numberOfPossibleIndices =
      range[1] - range[0] - previousIndices.size();
  static std::random_device generator;
  std::uniform_int_distribution<std::size_t> distribution(
      0, numberOfPossibleIndices - 1);

  std::size_t ind = distribution(generator); // New random position

  // Turn the random position into a previously not used index in the range.

  std::size_t newIndex = range[0];

  std::size_t count = 0;
  while (count <= ind) {
    if (!previousIndices.count(newIndex)) {
      count++;
      if (count <= ind)
        newIndex++;
      continue;
    }
    newIndex++;
  }
  previousIndices.insert(newIndex);
  return newIndex;
}
}

#endif
