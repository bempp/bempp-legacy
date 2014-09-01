// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_IMPL_HPP
#define HMAT_HMATRIX_IMPL_HPP

#include "hmatrix.hpp"
#include "hmatrix_data.hpp"
#include "hmatrix_dense_data.hpp"

#include <algorithm>

namespace hmat {

template <typename ValueType, int N>
HMatrix<ValueType, N>::HMatrix(
    const shared_ptr<BlockClusterTree<N>> &blockClusterTree)
    : m_blockClusterTree(blockClusterTree) {}

template <typename ValueType, int N>
HMatrix<ValueType, N>::HMatrix(
    const shared_ptr<BlockClusterTree<N>> &blockClusterTree,
    const HMatrixCompressor<ValueType, N> &hMatrixCompressor)
    : HMatrix<ValueType, N>(blockClusterTree) {
  initialize(hMatrixCompressor);
}

template <typename ValueType, int N>
std::size_t HMatrix<ValueType, N>::rows() const {
  return m_blockClusterTree->rows();
}

template <typename ValueType, int N>
std::size_t HMatrix<ValueType, N>::columns() const {
  return m_blockClusterTree->columns();
}

template <typename ValueType, int N>
void HMatrix<ValueType, N>::initialize(
    const HMatrixCompressor<ValueType, N> &hMatrixCompressor) {

  reset();

  for (auto node : m_blockClusterTree->leafNodes()) {
    shared_ptr<HMatrixData<ValueType>> nodeData;
    hMatrixCompressor.compressBlock(*node, nodeData);
    m_hMatrixData[node] = nodeData;
  }
}
template <typename ValueType, int N> void HMatrix<ValueType, N>::reset() {
  m_hMatrixData.clear();
}

template <typename ValueType, int N>
bool HMatrix<ValueType, N>::isInitialized() const {
  return (!m_hMatrixData.empty());
}

template <typename ValueType, int N>
arma::Mat<ValueType>
HMatrix<ValueType, N>::permuteMatToHMatDofs(const arma::Mat<ValueType> &mat,
                                            RowColSelector rowOrColumn) const {

  arma::Mat<ValueType> permutedDofs(mat.n_rows, mat.n_cols);

  shared_ptr<const ClusterTree<N>> clusterTree;

  if (rowOrColumn == ROW)
    clusterTree = m_blockClusterTree->rowClusterTree();
  else
    clusterTree = m_blockClusterTree->columnClusterTree();

  if (clusterTree->numberOfDofs() != mat.n_rows)
    throw std::runtime_error("HMatrix::permuteMatToHMatDofs: "
                             "Input matrix has wrong number of rows.");

  for (std::size_t i = 0; i < mat.n_rows; ++i) {
    auto permutedIndex = clusterTree->mapOriginalDofToHMatDof(i);
    for (std::size_t j = 0; j < mat.n_cols; ++j)
      permutedDofs(permutedIndex, j) = mat(i, j);
  }

  return permutedDofs;
}

template <typename ValueType, int N>
arma::Mat<ValueType> HMatrix<ValueType, N>::permuteMatToOriginalDofs(
    const arma::Mat<ValueType> &mat, RowColSelector rowOrColumn) const {

  arma::Mat<ValueType> originalDofs(mat.n_rows, mat.n_cols);

  shared_ptr<const ClusterTree<N>> clusterTree;

  if (rowOrColumn == ROW)
    clusterTree = m_blockClusterTree->rowClusterTree();
  else
    clusterTree = m_blockClusterTree->columnClusterTree();

  if (clusterTree->numberOfDofs() != mat.n_rows)
    throw std::runtime_error("HMatrix::permuteMatToOriginalDofs: "
                             "Input matrix has wrong number of rows.");

  for (std::size_t i = 0; i < mat.n_rows; ++i) {
    auto originalIndex = clusterTree->mapHMatDofToOriginalDof(i);
    for (std::size_t j = 0; j < mat.n_cols; ++j)
      originalDofs(originalIndex, j) = mat(i, j);
  }

  return originalDofs;
}

template <typename ValueType, int N>
void HMatrix<ValueType, N>::apply(const arma::Mat<ValueType> &X,
                                  arma::Mat<ValueType> &Y, TransposeMode trans,
                                  ValueType alpha, ValueType beta) const {

  if (beta == ValueType(0))
    Y.zeros();
  else
    Y *= beta;

  arma::Mat<ValueType> xPermuted;
  arma::Mat<ValueType> yPermuted;

  if (trans == TransposeMode::NOTRANS) {

    xPermuted = permuteMatToHMatDofs(X, COL);
    yPermuted = permuteMatToHMatDofs(Y, ROW);
  } else {
    xPermuted = permuteMatToHMatDofs(X, ROW);
    yPermuted = permuteMatToHMatDofs(Y, COL);
  }

  std::for_each(begin(m_hMatrixData), end(m_hMatrixData),
                [trans, alpha, beta, &xPermuted, &yPermuted, this](
                    const std::pair<shared_ptr<BlockClusterTreeNode<N>>,
                                    shared_ptr<HMatrixData<ValueType>>>
                        elem) {

    IndexRangeType inputRange;
    IndexRangeType outputRange;
    if (trans == TransposeMode::NOTRANS) {
      inputRange = elem.first->data().columnClusterTreeNode->data().indexRange;
      outputRange = elem.first->data().rowClusterTreeNode->data().indexRange;
    } else {
      inputRange = elem.first->data().rowClusterTreeNode->data().indexRange;
      outputRange = elem.first->data().columnClusterTreeNode->data().indexRange;
    }

    arma::subview<ValueType> xData =
        xPermuted.rows(inputRange[0], inputRange[1] - 1);
    arma::subview<ValueType> yData =
        yPermuted.rows(outputRange[0], outputRange[1] - 1);
    elem.second->apply(xData, yData, trans, alpha, 1);
  });

  Y = this->permuteMatToOriginalDofs(yPermuted, ROW);
}
}

#endif
