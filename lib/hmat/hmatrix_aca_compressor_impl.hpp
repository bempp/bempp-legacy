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

  A.resize(numberOfRows,0);
  B.resize(0,numberOfColumns);

  Vector<ValueType> randCompareVec = Vector<ValueType>::Random(B.cols());
  randCompareVec /= randCompareVec.norm();
  Vector<ValueType> randCompareResVec = Vector<ValueType>::Zero(A.rows());

  std::set<std::size_t> previousRowIndices;
  std::set<std::size_t> previousColumnIndices;

  std::size_t iterationLimit =
      std::min(static_cast<std::size_t>(m_maxRank),
               std::min(numberOfRows, numberOfColumns));

  std::size_t rankCount = 0;
  std::ptrdiff_t maxRowInd;
  std::ptrdiff_t maxColInd;

  IndexRangeType rowIndexRange;
  IndexRangeType columnIndexRange;


  for (int i = 0; i < iterationLimit; ++i) {

    std::size_t rowIndex = randomIndex(rowClusterRange, previousRowIndices);

    // Compute complete row

    Matrix<ValueType> newRow;
    Matrix<ValueType> newCol;

    rowIndexRange = {{rowIndex, rowIndex + 1}};
    columnIndexRange = columnClusterRange;
    
    m_dataAccessor.computeMatrixBlock(rowIndexRange, columnIndexRange,
                                      blockClusterTreeNode, newRow);

    if (A.cols()>0)
      newRow -= A.row((rowIndex-rowClusterRange[0]))*B;

    auto val = newRow.cwiseAbs().maxCoeff(&maxRowInd, &maxColInd);

    if (val < 1E-14)
      continue;

    auto pivot = newRow(0, maxColInd);

    newRow /= pivot;

    // Now evaluate column

    maxColInd += columnClusterRange[0]; // Map back to original variables

    rowIndexRange = rowClusterRange;
    columnIndexRange = {{boost::numeric_cast<std::size_t>(maxColInd),
                         boost::numeric_cast<std::size_t>(maxColInd) + 1}};

    m_dataAccessor.computeMatrixBlock(rowIndexRange, columnIndexRange,
       blockClusterTreeNode, newCol);

    if (B.rows()>0)
      newCol -= A*B.col((maxColInd-columnClusterRange[0]));
    
    // Enlarge A and B
    Matrix<ValueType> Atmp(A.rows(),A.cols()+1);
    Matrix<ValueType> Btmp(B.rows()+1,B.cols());
    Atmp.leftCols(A.cols()) = A;
    Btmp.topRows(B.rows()) = B;
    A.swap(Atmp);
    B.swap(Btmp);

    A.rightCols(1) = newCol;
    B.bottomRows(1) = newRow;

    //auto frobeniusNorm = hMatrixData->frobeniusNorm();

    randCompareResVec += A.rightCols(1)*(B.bottomRows(1)*randCompareVec);
    auto normEst = randCompareResVec.norm();

    if (newCol.norm() * newRow.norm() < m_eps * normEst){
      //std::cout << rankCount << std::endl;
      break;
    }
    

  }
  // if (isnan(A) || isnan(B))
  //  throw std::runtime_error("NaN detected before end of compress.");

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
