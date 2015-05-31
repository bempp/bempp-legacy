// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_IMPL_HPP
#define HMAT_HMATRIX_IMPL_HPP

#include "hmatrix.hpp"
#include "hmatrix_data.hpp"
#include "hmatrix_dense_data.hpp"
#include <tbb/parallel_for_each.h>
#include <tbb/task_group.h>

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

  typedef decltype(m_blockClusterTree->root()) node_t;
  
  std::function<void(const node_t& node)> compressFun = 
   [&](const node_t& node){
    if (node->isLeaf())
    {
     shared_ptr<HMatrixData<ValueType>> nodeData;
     hMatrixCompressor.compressBlock(*node, nodeData);
     m_hMatrixData[node] = nodeData;
    }
    else 
    {
      tbb::task_group g;
      g.run([&]{compressFun(node->child(0));});
      g.run([&]{compressFun(node->child(1));});
      g.run([&]{compressFun(node->child(2));});
      g.run([&]{compressFun(node->child(3));});
      g.wait();
    }
  };
            

  compressFun(m_blockClusterTree->root());


  // The following can be removed for C++14 compilers that
  // support polymorphic lambda expressions.
  //typedef decltype(*begin(leafNodes)) node_t;

  //tbb::parallel_for_each(begin(leafNodes),end(leafNodes),[&](const node_t& node){
  //  shared_ptr<HMatrixData<ValueType>> nodeData;
  //  hMatrixCompressor.compressBlock(*node, nodeData);
  //  m_hMatrixData[node] = nodeData;
  //});
}

template <typename ValueType, int N> void HMatrix<ValueType, N>::reset() {
  m_hMatrixData.clear();
}

template <typename ValueType, int N>
bool HMatrix<ValueType, N>::isInitialized() const {
  return (!m_hMatrixData.empty());
}

template <typename ValueType, int N>
shared_ptr<const BlockClusterTree<N>> 
HMatrix<ValueType, N>::blockClusterTree() const {
  return this->m_blockClusterTree;
}

template <typename ValueType, int N>
shared_ptr<const HMatrixData<ValueType>> 
HMatrix<ValueType, N>::data(
    shared_ptr<const BlockClusterTreeNode<N>> node) const {
  return this->m_hMatrixData.at(const_pointer_cast<BlockClusterTreeNode<N>>(
      node));
}

template <typename ValueType, int N>
Matrix<ValueType>
HMatrix<ValueType, N>::permuteMatToHMatDofs(const Matrix<ValueType> &mat,
                                            RowColSelector rowOrColumn) const {

  Matrix<ValueType> permutedDofs(mat.rows(), mat.cols());

  shared_ptr<const ClusterTree<N>> clusterTree;

  if (rowOrColumn == ROW)
    clusterTree = m_blockClusterTree->rowClusterTree();
  else
    clusterTree = m_blockClusterTree->columnClusterTree();

  if (clusterTree->numberOfDofs() != mat.rows())
    throw std::runtime_error("HMatrix::permuteMatToHMatDofs: "
                             "Input matrix has wrong number of rows.");

  for (std::size_t i = 0; i < mat.rows(); ++i) {
    auto permutedIndex = clusterTree->mapOriginalDofToHMatDof(i);
    for (std::size_t j = 0; j < mat.cols(); ++j)
      permutedDofs(permutedIndex, j) = mat(i, j);
  }

  return permutedDofs;
}

template <typename ValueType, int N>
Matrix<ValueType> HMatrix<ValueType, N>::permuteMatToOriginalDofs(
    const Matrix<ValueType> &mat, RowColSelector rowOrColumn) const {

  Matrix<ValueType> originalDofs(mat.rows(), mat.cols());

  shared_ptr<const ClusterTree<N>> clusterTree;

  if (rowOrColumn == ROW)
    clusterTree = m_blockClusterTree->rowClusterTree();
  else
    clusterTree = m_blockClusterTree->columnClusterTree();

  if (clusterTree->numberOfDofs() != mat.rows())
    throw std::runtime_error("HMatrix::permuteMatToOriginalDofs: "
                             "Input matrix has wrong number of rows.");

  for (std::size_t i = 0; i < mat.rows(); ++i) {
    auto originalIndex = clusterTree->mapHMatDofToOriginalDof(i);
    for (std::size_t j = 0; j < mat.cols(); ++j)
      originalDofs(originalIndex, j) = mat(i, j);
  }

  return originalDofs;
}

template <typename ValueType, int N>
void HMatrix<ValueType, N>::apply(const Eigen::Ref<Matrix<ValueType>>& X,
                                  Eigen::Ref<Matrix<ValueType>> Y, TransposeMode trans,
                                  ValueType alpha, ValueType beta) const {

  if (beta == ValueType(0))
    Y.setZero();
  else
    Y *= beta;

  typedef decltype(m_blockClusterTree->root()) node_t;

  Matrix<ValueType> xPermuted;
  Matrix<ValueType> yPermuted;

  if (trans == TransposeMode::NOTRANS) {

    xPermuted = permuteMatToHMatDofs(X, COL);
    yPermuted = permuteMatToHMatDofs(Y, ROW);
  } else {
    xPermuted = permuteMatToHMatDofs(X, ROW);
    yPermuted = permuteMatToHMatDofs(Y, COL);
  }


  std::function<void(const node_t&, const Eigen::Ref<Matrix<ValueType>>&,
      Eigen::Ref<Matrix<ValueType>>)> applyFun = [&](
        const node_t& node, const Eigen::Ref<Matrix<ValueType>>& x_in,
        Eigen::Ref<Matrix<ValueType>> y_inout)
      {
        if (node->isLeaf())
        {
          m_hMatrixData[node]->apply(x_in,y_inout,trans,alpha,1.0);
        }
        else
        {



      };

    


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

    Matrix<ValueType> xData = xPermuted.block(inputRange[0],0,inputRange[1]-inputRange[0],xPermuted.cols());
    Matrix<ValueType> yData = yPermuted.block(outputRange[0],0,outputRange[1]-outputRange[0],yPermuted.cols());


    elem.second->apply(xData, yData, trans, alpha, 1);

    yPermuted.block(outputRange[0],0,outputRange[1]-outputRange[0],yData.cols()) = yData;
  });

  Y = this->permuteMatToOriginalDofs(yPermuted, ROW);
}
}

#endif
