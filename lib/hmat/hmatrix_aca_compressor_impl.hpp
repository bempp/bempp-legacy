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
bool HMatrixAcaCompressor<ValueType, N>::selectMinPivot(
    const Matrix<ValueType> &vec,
    const std::vector<std::size_t> &approximationCount, std::size_t &pivot,
    ModeType modus) const {

  int n = std::max(vec.rows(), vec.cols());
  pivot = n;
  int row = 0;
  int col = 0;
  double minVal = std::numeric_limits<double>::max();
  for (int i = 0; i < n; ++i) {
    if (modus == ModeType::ROW)
      col = i;
    else
      row = i;
    if (approximationCount[i] == 0 && std::abs(vec(row, col)) < minVal) {
      minVal = std::abs(vec(row, col));
      pivot = i;
    }
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
    std::vector<std::size_t> &colApproxCounter, ModeType mode,
    double zeroTol) const {

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
    // is associated with zero for the search.
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
    // is associated with zero for the search.
    ValueType tmp = 0;
    std::swap(tmp, row(0, columnIndex - columnClusterRange[0]));
    row.cwiseAbs().maxCoeff(&maxRowInd, &maxColInd);
    std::swap(tmp, row(0, columnIndex - columnClusterRange[0]));
    nextPivot = boost::numeric_cast<std::size_t>(maxColInd);
  }

  rowApproxCounter[rowIndex - rowClusterRange[0]] += 1;
  colApproxCounter[columnIndex - columnClusterRange[0]] += 1;

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

  Matrix<ValueType> origRow, origCol;

  getBlockClusterTreeNodeDimensions(blockClusterTreeNode, rowClusterRange,
                                    columnClusterRange, numberOfRows,
                                    numberOfColumns);

  hMatrixData.reset(new HMatrixLowRankData<ValueType>());

  A.resize(numberOfRows, 0);
  B.resize(0, numberOfColumns);

  std::vector<size_t> rowApproxCounter(numberOfRows);
  std::vector<size_t> colApproxCounter(numberOfColumns);

  double blockNorm = 0;
  std::size_t maxIterations = std::min(static_cast<std::size_t>(m_maxRank),
                                       std::min(numberOfRows, numberOfColumns));

  bool finished = false;
  AcaAlgorithmStateType state = AcaAlgorithmStateType::START;
  AcaStatusType acaStatus;
  std::size_t nextPivot = 0;
  ModeType mode = ModeType::ROW;

  bool rowTrialPassed = false;
  bool columnTrialPassed = false;

  // Count how often row or column trial has been executed
  std::size_t rowTrialCount = 0;
  std::size_t columnTrialCount = 0;

  // constants for maximum row and column trial counts

  const std::size_t MAX_ROW_TRIAL_COUNT = 1;
  const std::size_t MAX_COLUMN_TRIAL_COUNT = 1;

  bool rowModeZeroTermination = false;
  bool columnModeZeroTermination = false;

  while (!finished) {
//  if(!finished) {
    // First run the ACA
    acaStatus = aca(blockClusterTreeNode, nextPivot, A, B, maxIterations,
                    rowApproxCounter, colApproxCounter, origRow, origCol,
                    blockNorm, m_eps, 1E-15, mode);
    // Now test the status for different cases
    if (acaStatus == AcaStatusType::RANK_LIMIT_REACHED) {
      finished = true; // Approximation does not work. So just stop.
    }
    // The next should only happen if we are in ROW_TRIAL or COLUMN_TRIAL mode
    if (acaStatus == AcaStatusType::CONVERGED_WITHOUT_ITERATION) {
      if (state == AcaAlgorithmStateType::START)
        state = AcaAlgorithmStateType::ROW_TRIAL;
      else if (state == AcaAlgorithmStateType::ROW_TRIAL) {
        rowTrialPassed = true;
        state = AcaAlgorithmStateType::COLUMN_TRIAL;
      } else if (state == AcaAlgorithmStateType::COLUMN_TRIAL) {
        columnTrialPassed = true;
        state = AcaAlgorithmStateType::ROW_TRIAL;
      }
      // If both, the row trial and column trial are passed, finish.
      if (rowTrialPassed && columnTrialPassed)
        finished = true;
    }
    // The next is the normal case after running ACA
    if (acaStatus == AcaStatusType::CONVERGED_AFTER_ITERATION ||
        acaStatus == AcaStatusType::ZERO_TERMINATION_AFTER_ITERATION) {
      // Have to rexecute row trial and column trial
      rowTrialPassed = false;
      columnTrialPassed = false;
      rowModeZeroTermination = false;
      columnModeZeroTermination = false;
      // Choose which trial to go into
      state = (mode == ModeType::ROW) ? AcaAlgorithmStateType::ROW_TRIAL
                                      : AcaAlgorithmStateType::COLUMN_TRIAL;
    }
    if (acaStatus == AcaStatusType::ZERO_TERMINATION_WITHOUT_ITERATION) {
      if (mode == ModeType::ROW) {
        rowModeZeroTermination = true;
        state = AcaAlgorithmStateType::ROW_TRIAL;
      } else {
        columnModeZeroTermination = true;
        state = AcaAlgorithmStateType::COLUMN_TRIAL;
      }
      if (rowModeZeroTermination && columnModeZeroTermination)
        finished = true;
    }
    if (!finished && state == AcaAlgorithmStateType::ROW_TRIAL) {
      if ((rowTrialCount < MAX_ROW_TRIAL_COUNT) &&
          selectMinPivot(origRow, colApproxCounter, nextPivot, ModeType::ROW)) {
        mode = ModeType::COL;
        rowTrialCount++;
      } else {
        // All columns have been approximated
        finished = true;
      }
    }
    if (!finished && state == AcaAlgorithmStateType::COLUMN_TRIAL) {
      if ((columnTrialCount < MAX_COLUMN_TRIAL_COUNT) &&
          selectMinPivot(origCol, rowApproxCounter, nextPivot, ModeType::COL)) {
        mode = ModeType::ROW;
        columnTrialCount++;
      } else {
        // All columns have been approximated
        finished = true;
      }
    }
  }
  static_cast<HMatrixLowRankData<ValueType> *>(hMatrixData.get())->A().swap(A);
  static_cast<HMatrixLowRankData<ValueType> *>(hMatrixData.get())->B().swap(B);
}

template <typename ValueType, int N>
typename HMatrixAcaCompressor<ValueType, N>::AcaStatusType
HMatrixAcaCompressor<ValueType, N>::aca(
    const BlockClusterTreeNode<N> &blockClusterTreeNode, std::size_t startPivot,
    Matrix<ValueType> &A, Matrix<ValueType> &B, std::size_t &maxIterations,
    std::vector<size_t> &rowApproxCounter,
    std::vector<size_t> &colApproxCounter, Matrix<ValueType> &origRow,
    Matrix<ValueType> &origCol, double &blockNorm, double eps, double zeroTol,
    ModeType mode) const {

  IndexRangeType rowClusterRange;
  IndexRangeType columnClusterRange;
  std::size_t numberOfRows;
  std::size_t numberOfColumns;

  getBlockClusterTreeNodeDimensions(blockClusterTreeNode, rowClusterRange,
                                    columnClusterRange, numberOfRows,
                                    numberOfColumns);

  bool converged = false;
  std::size_t nextPivot = startPivot;
  std::size_t iterationCount = 0;

  CrossStatusType crossStatus;
  Matrix<ValueType> row, col;

  while (maxIterations > 0) {
//  if (maxIterations > 0) {
    crossStatus = computeCross(blockClusterTreeNode, A, B, nextPivot, origRow,
                               origCol, row, col, rowApproxCounter,
                               colApproxCounter, mode, zeroTol);
    if (crossStatus == CrossStatusType::ZERO)
      return (iterationCount == 0)
                 ? AcaStatusType::ZERO_TERMINATION_WITHOUT_ITERATION
                 : AcaStatusType::CONVERGED_AFTER_ITERATION;
    converged = checkConvergence(row, col, eps * blockNorm);
    if (converged)
      break;
    else {
      blockNorm = updateLowRankBlocksAndNorm(row, col, A, B, blockNorm);
      maxIterations--; // Putting it here implies that cross computation for
                       // convergence testing does not count towards the max
                       // number of iterations.
    }
    iterationCount++;
  }

  if (converged)
    return (iterationCount == 0) ? AcaStatusType::CONVERGED_WITHOUT_ITERATION
                                 : AcaStatusType::CONVERGED_AFTER_ITERATION;
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
