// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_ACA_COMPRESSOR_IMPL_HPP
#define HMAT_HMATRIX_ACA_COMPRESSOR_IMPL_HPP

#include "hmatrix_aca_compressor.hpp"
#include <boost/numeric/conversion/cast.hpp>
#include "hmatrix_low_rank_data.hpp"
#include "eigen_fwd.hpp"
#include "scalar_traits.hpp"
#include <random>
#include <complex>
#include <cmath>
#include <algorithm>
#include <limits>

namespace hmat {

template <typename ValueType, int N>
bool HMatrixAcaCompressor<ValueType, N>::selectPivotWithMinNorm2(
    const std::vector<double> &norms,
    const std::vector<int> &approximationCount, int &pivot) {

  int n = norms.size();
  pivot = n;
  double minNorm = std::numeric_limits<double>::max();
  for (int i = 0; i < n; ++i)
    if (approximationCount[i] == 0 && norms[i] < minNorm) {
      minNorm = norms[i];
      pivot = i;
    }
  return (pivot != n);
}

template <typename ValueType, int N>
typename HMatrixAcaCompressor<ValueType, N>::CrossStatusType
HMatrixAcaCompressor<ValueType, N>::computeCross(
    const BlockClusterTreeNode<N> &blockClusterTreeNode,
    const Matrix<ValueType> &A, const Matrix<ValueType> &B,
    std::size_t &nextPivot, Matrix<ValueType> &origRow,
    Matrix<ValueType> &origCol, Matrix<ValueType> &row, Matrix<ValueType> &col,
    std::vector<std::size_t> &rowApproxCounter,
    std::vector<std::size_t> &colApproxCounter, std::vector<double> &rowNorms,
    std::vector<double> &colNorms, ModeType mode, double zeroTol) const {

  const auto &rowClusterRange =
      blockClusterTreeNode.data().rowClusterTreeNode->data().indexRange;

  const auto &columnClusterRange =
      blockClusterTreeNode.data().columnClusterTreeNode->data().indexRange;

  std::size_t rowIndex;
  std::size_t columnIndex;

  std::ptrdiff_t maxRowInd;
  std::ptrdiff_t maxColInd;

  IndexRangeType rowIndexRange;
  IndexRangeType columnIndexRange;

  ValueType pivotValue;

  if (mode == ModeType::ROW) {

    rowIndex = rowClusterRange[0] + nextPivot;
    rowIndexRange = {{rowIndex, rowIndex + 1}};
    m_dataAccessor.computeMatrixBlock(rowIndexRange, columnClusterRange,
                                      blockClusterTreeNode, origRow);
    row = origRow;
    if (A.cols() > 0 && B.rows() > 0)
      row -= A.row(rowIndex - rowClusterRange[0]) * B;

    if (row.norm() <= zeroTol)
      return CrossStatusType::ZERO;

    row.cwiseAbs().maxCoeff(&maxRowInd, &maxColInd);
    pivotValue = row(0, maxColInd);
    columnIndex =
        boost::numeric_cast<std::size_t>(maxColInd) + columnClusterRange[0];
    columnIndexRange = {{columnIndex, columnIndex + 1}};
    m_dataAccessor.computeMatrixBlock(rowClusterRange, columnIndexRange,
                                      blockClusterTreeNode, origCol);

    col = origCol;
    if (A.cols() > 0 && B.rows() > 0)
      col -= A * B.col(columnIndex - columnClusterRange[0]);

    row /= pivotValue;

    // Find next pivot, make sure that current rowIndex
    // won't be taken.
    ValueType tmp = 0;
    std::swap(tmp, col(rowIndex - rowClusterRange[0], 0));
    col.cwiseAbs().maxCoeff(&maxRowInd, &maxColInd);
    std::swap(tmp, col(rowIndex - rowClusterRange[0], 0));
    nextPivot = boost::numeric_cast<std::size_t>(maxRowInd);

  }

  else {
    columnIndex = columnClusterRange[0] + nextPivot;
    columnIndexRange = {{columnIndex, columnIndex + 1}};
    m_dataAccessor.computeMatrixBlock(rowClusterRange, columnIndexRange,
                                      blockClusterTreeNode, origCol);
    col = origCol;
    if (A.cols() > 0 && B.rows() > 0)
      col -= A * B.col(columnIndex - columnClusterRange[0]);

    if (col.norm() <= zeroTol)
      return CrossStatusType::ZERO;

    col.cwiseAbs().maxCoeff(&maxRowInd, &maxColInd);
    pivotValue = col(maxRowInd, 0);
    rowIndex = boost::numeric_cast<std::size_t>(maxRowInd) + rowClusterRange[0];
    rowIndexRange = {{rowIndex, rowIndex + 1}};
    m_dataAccessor.computeMatrixBlock(rowIndexRange, columnClusterRange,
                                      blockClusterTreeNode, origRow);

    row = origRow;
    if (A.cols() > 0 && B.rows() > 0)
      row -= A.row(rowIndex - rowClusterRange[0]) * B;

    col /= pivotValue;

    // Find next pivot, make sure that current columnIndex
    // won't be taken.
    ValueType tmp = 0;
    std::swap(tmp, row(0, columnIndex - columnClusterRange[0]));
    row.cwiseAbs().maxCoeff(&maxRowInd, &maxColInd);
    std::swap(tmp, row(0, columnIndex - columnClusterRange[0]));
    nextPivot = boost::numeric_cast<std::size_t>(maxColInd);
  }

  rowApproxCounter[rowIndex - rowClusterRange[0]] += 1;
  colApproxCounter[columnIndex - columnClusterRange[0]] += 1;
  rowNorms[rowIndex - rowClusterRange[0]] = origRow.norm();
  colNorms[columnIndex - columnClusterRange[0]] = origCol.norm();

  return CrossStatusType::SUCCESS;
}

template <typename ValueType, int N>
double HMatrixAcaCompressor<ValueType, N>::updateLowRankBlocksAndNorm(
    const Matrix<ValueType> &newRow, const Matrix<ValueType> &newCol,
    Matrix<ValueType> &A, Matrix<ValueType> &B, double currentBlockNorm) const {

  // Compute norm update

  double newNormSquared = newRow.squaredNorm() * newCol.squaredNorm();

  if (A.cols() > 0 && B.rows() > 0) {
    double val1 =
        2 * (newRow * (B.adjoint() * (A.adjoint() * newCol))).real()(0, 0);
    double val2 = currentBlockNorm * currentBlockNorm;
    newNormSquared += val1 + val2;
  }

  // Now increase size of blocks
  Matrix<ValueType> Atmp(A.rows(), A.cols() + 1);
  Matrix<ValueType> Btmp(B.rows() + 1, B.cols());
  Atmp.leftCols(A.cols()) = A;
  Btmp.topRows(B.rows()) = B;
  A.swap(Atmp);
  B.swap(Btmp);

  A.rightCols(1) = newCol;
  B.bottomRows(1) = newRow;

  // Return 0 if due to rounding errors newNormSquared is smaller 0.
  return std::sqrt((newNormSquared > 0) ? newNormSquared : 0);
}

template <typename ValueType, int N>
bool HMatrixAcaCompressor<ValueType, N>::checkConvergence(
    const Matrix<ValueType> &row, const Matrix<ValueType> &col,
    double tol) const {

  return row.norm() * col.norm() < tol;
}

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

  Matrix<ValueType> A;
  Matrix<ValueType> B;

  getBlockClusterTreeNodeDimensions(blockClusterTreeNode, rowClusterRange,
                                    columnClusterRange, numberOfRows,
                                    numberOfColumns);

  hMatrixData.reset(new HMatrixLowRankData<ValueType>());

  A.resize(numberOfRows, 0);
  B.resize(0, numberOfColumns);

  std::vector<size_t> rowApproxCounter(numberOfRows);
  std::vector<size_t> colApproxCounter(numberOfColumns);

  std::vector<double> rowNorms(numberOfRows);
  std::vector<double> colNorms(numberOfColumns);

  double blockNorm = 0;
  std::size_t iterationLimit =
      std::min(static_cast<std::size_t>(m_maxRank),
               std::min(numberOfRows, numberOfColumns));

  AcaStatusType status = aca(blockClusterTreeNode, 0, A, B, iterationLimit,
                             rowApproxCounter, colApproxCounter, rowNorms,
                             colNorms, blockNorm, m_eps, 1E-15, ModeType::ROW);

  static_cast<HMatrixLowRankData<ValueType> *>(hMatrixData.get())->A().swap(A);
  static_cast<HMatrixLowRankData<ValueType> *>(hMatrixData.get())->B().swap(B);
}

template <typename ValueType, int N>
typename HMatrixAcaCompressor<ValueType, N>::AcaStatusType
HMatrixAcaCompressor<ValueType, N>::aca(
    const BlockClusterTreeNode<N> &blockClusterTreeNode, std::size_t startPivot,
    Matrix<ValueType> &A, Matrix<ValueType> &B, std::size_t maxIterations,
    std::vector<size_t> &rowApproxCounter,
    std::vector<size_t> &colApproxCounter, std::vector<double> &rowNorms,
    std::vector<double> &colNorms, double &blockNorm, double eps,
    double zeroTol, ModeType mode) const {

  IndexRangeType rowClusterRange;
  IndexRangeType columnClusterRange;
  std::size_t numberOfRows;
  std::size_t numberOfColumns;

  getBlockClusterTreeNodeDimensions(blockClusterTreeNode, rowClusterRange,
                                    columnClusterRange, numberOfRows,
                                    numberOfColumns);

  bool converged = false;
  std::size_t nextPivot = startPivot;
  std::size_t rankCount = 0;

  CrossStatusType crossStatus;
  Matrix<ValueType> origRow, origCol, row, col;

  while (!converged && rankCount < maxIterations) {
    crossStatus = computeCross(
        blockClusterTreeNode, A, B, nextPivot, origRow, origCol, row, col,
        rowApproxCounter, colApproxCounter, rowNorms, colNorms, mode, zeroTol);
    if (crossStatus == CrossStatusType::ZERO)
      return AcaStatusType::ZERO_TERMINATION;

    blockNorm = updateLowRankBlocksAndNorm(row, col, A, B, blockNorm);
    converged = checkConvergence(row, col, eps * blockNorm);
    rankCount++;
  }

  if (converged)
    return AcaStatusType::SUCCESS;
  else
    return AcaStatusType::RANK_LIMIT_REACHED;
}

template <typename ValueType, int N>
HMatrixAcaCompressor<ValueType, N>::HMatrixAcaCompressor(
    const DataAccessor<ValueType, N> &dataAccessor, double eps,
    unsigned int maxRank)
    : m_dataAccessor(dataAccessor), m_eps(eps), m_maxRank(maxRank),
      m_hMatrixDenseCompressor(dataAccessor) {}

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
