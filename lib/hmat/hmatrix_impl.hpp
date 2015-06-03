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
    : m_blockClusterTree(blockClusterTree),
      m_numberOfDenseBlocks(0),
      m_numberOfLowRankBlocks(0),
      m_memSizeKb(0.0) {}

template <typename ValueType, int N>
HMatrix<ValueType, N>::HMatrix(
    const shared_ptr<BlockClusterTree<N>> &blockClusterTree,
    const HMatrixCompressor<ValueType, N> &hMatrixCompressor,
    bool coarsening, double coarsening_accuracy)
    : HMatrix<ValueType, N>(blockClusterTree) {
  initialize(hMatrixCompressor, coarsening, coarsening_accuracy);
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
double HMatrix<ValueType, N>::frobeniusNorm() const {
  return frobeniusNorm_impl(m_blockClusterTree->root());
}

template <typename ValueType, int N>
int HMatrix<ValueType, N>::numberOfDenseBlocks() const {
  return m_numberOfDenseBlocks;
}

template <typename ValueType, int N>
int HMatrix<ValueType, N>::numberOfLowRankBlocks() const {
  return m_numberOfLowRankBlocks;
}

template <typename ValueType, int N>
int HMatrix<ValueType, N>::numberOfBlocks() const {
  return m_numberOfLowRankBlocks+m_numberOfDenseBlocks;
}

template <typename ValueType, int N>
double HMatrix<ValueType, N>::memSizeKb() const {
  return m_memSizeKb;
}

template <typename ValueType, int N>
void HMatrix<ValueType, N>::initialize(
    const HMatrixCompressor<ValueType, N> &hMatrixCompressor,
    bool coarsening, double coarsening_accuracy) {

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
      g.run_and_wait([&]{compressFun(node->child(3));});
    }


    // Check if coarsening is possible
    
    bool children_are_leafs = true;
    bool children_are_admissible = true;

    for (int i = 0; i < 4; ++i){
      children_are_leafs = node->child(i)->isLeaf() && children_are_leafs;
      children_are_admissible = node->child(i)->admissible() && children_are_admissible;
    }
    
    bool start_coarsen = coarsening && children_are_leafs && children_are_admissible;

    if (start_coarsen){

      // get max rank of children
      
      int rank = 0;
      for (int i =0; i < 4; ++i)
        rank = std::max(rank,node->child(i)->data().rank());

      int p = 4;
      int k = rank+p;

      bool accept = false;

      IndexRangeType rowIndexRange = node->data().rowClusterTree->data().indexRange;
      IndexRangeType colIndexRange = node->data().columnClusterTree->data().indexRange;
      int rows = rowIndexRange[1]-rowIndexRange[0];
      int cols = colIndexrange[1]-colIndexRange[0];

      double norm = frobeniusNorm_impl(node);

      int count = 0;
      while (!accept) {
        if (k>cols || count > 1) break;

        Matrix<ValueType> Z = Matrix<ValueType>::Random(cols,k);
        Matrix<ValueType> V = Matrix<ValueType>::Zero(rows,k);

        apply_impl(node,Eigen::Ref<Matrix<ValueType>>(Z),Eigen::Ref<Matrix<ValueType>>(V),
            NOTRANS);
        Eigen::HouseholderQR qr(V);
        Matrix<ValueType> Q = qr.householderQ();
        Matrix<ValueType> B(cols,k);
        apply_impl(node,Eigen::Ref<Matrix<ValueType>>(Q),Eigen::Ref<Matrix<ValueType>>(B),
            CONJTRANS);
        Eigen::JacobiSVD svd(B);
        auto singularValues = svd.singularValues();
        singularValues *= 1./norm;
        int index;
        for (index = 0; index != singularValues.rows(); ++index)
          if (singularValues(index)<coarsening_accuracy) break;

        if (index!=singularValues.rows()){
          // Create node and get out of loop
          

          

          
        }



        



        

        



      } 

      




    }

  };
    

  compressFun(m_blockClusterTree->root());

  for (auto& elem: m_hMatrixData){
     if (!elem.second) continue;
     if (elem.second->type()==DENSE)
       m_numberOfDenseBlocks++;
     else
       m_numberOfLowRankBlocks++;
     m_memSizeKb+=elem.second->memSizeKb();
  }

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
HMatrix<ValueType, N>::permuteMatToHMatDofs(const Eigen::Ref<Matrix<ValueType>> &mat,
                                            RowColSelector rowOrColumn) const {

  Matrix<ValueType> permutedDofs(mat.rows(), mat.cols());

  shared_ptr<const ClusterTree<N>> clusterTree;

  if (rowOrColumn == ROW)
    clusterTree = m_blockClusterTree->rowClusterTree();
  else
    clusterTree = m_blockClusterTree->columnClusterTree();

  for (std::size_t i = 0; i < mat.rows(); ++i) {
    auto permutedIndex = clusterTree->mapOriginalDofToHMatDof(i);
    for (std::size_t j = 0; j < mat.cols(); ++j)
      permutedDofs(permutedIndex, j) = mat(i, j);
  }

  return permutedDofs;
}

template <typename ValueType, int N>
Matrix<ValueType> HMatrix<ValueType, N>::permuteMatToOriginalDofs(
    const Eigen::Ref<Matrix<ValueType>> &mat, RowColSelector rowOrColumn) const {

  Matrix<ValueType> originalDofs(mat.rows(), mat.cols());

  shared_ptr<const ClusterTree<N>> clusterTree;

  if (rowOrColumn == ROW)
    clusterTree = m_blockClusterTree->rowClusterTree();
  else
    clusterTree = m_blockClusterTree->columnClusterTree();

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

  Matrix<ValueType> xPermuted;
  Matrix<ValueType> yPermuted;


  if (trans == TransposeMode::NOTRANS || trans == TransposeMode::CONJ) {

    xPermuted = alpha*permuteMatToHMatDofs(X, COL);
    yPermuted = permuteMatToHMatDofs(Y, ROW);
  } else {
    xPermuted = alpha*permuteMatToHMatDofs(X, ROW);
    yPermuted = permuteMatToHMatDofs(Y, COL);
  }

  apply_impl(this->m_blockClusterTree->root(),Eigen::Ref<Matrix<ValueType>>(xPermuted), 
      Eigen::Ref<Matrix<ValueType>>(yPermuted), trans);


  if (trans == TransposeMode::NOTRANS || trans == TransposeMode::CONJ)
    Y = this->permuteMatToOriginalDofs(yPermuted, ROW);
  else
    Y = this->permuteMatToOriginalDofs(yPermuted, COL);
  

    

  // Y = this->permuteMatToOriginalDofs(yPermuted, ROW);
}


template <typename ValueType, int N>
void HMatrix<ValueType, N>::apply_impl(const shared_ptr<BlockClusterTreeNode<N>>& node,
      const Eigen::Ref<Matrix<ValueType>> &X, Eigen::Ref<Matrix<ValueType>> Y,
      TransposeMode trans) const
{

  if (node->isLeaf()){
    this->m_hMatrixData.at(node)->apply(X,Y,trans,1,1);
  }
  else
  {

   auto child0= node->child(0);
   auto child1= node->child(1);
   auto child2= node->child(2);
   auto child3= node->child(3);

   IndexRangeType childRowRange = child0->data().rowClusterTreeNode->data().indexRange;
   IndexRangeType childColRange = child0->data().columnClusterTreeNode->data().indexRange;

   int rowSplit = childRowRange[1]-childRowRange[0];
   int colSplit = childColRange[1]-childColRange[0];

   IndexRangeType inputRange;
   IndexRangeType outputRange;

   int x_size = X.rows();
   int y_size = Y.rows();
   
    // Ugly hack since Eigen::Ref constructor does not like const objects
    Eigen::Ref<Matrix<ValueType>> x_no_const = const_cast<Eigen::Ref<Matrix<ValueType>>&>(X);

    int x_split;
    int y_split;

    if (trans == TransposeMode::NOTRANS || trans == TransposeMode::CONJ) 
    {
      x_split = colSplit;
      y_split = rowSplit;
    } 
    else 
    {
      x_split = rowSplit;
      y_split = colSplit;
    }

    Eigen::Ref<Matrix<ValueType>> xData1(x_no_const.topRows(x_split));
    Eigen::Ref<Matrix<ValueType>> xData2(x_no_const.bottomRows(X.rows()-x_split));

    Eigen::Ref<Matrix<ValueType>> yData1(Y.topRows(y_split));
    Eigen::Ref<Matrix<ValueType>> yData2(Y.bottomRows(Y.rows()-y_split));
    
    tbb::task_group g;

    if (trans == TransposeMode::NOTRANS || trans == TransposeMode::CONJ) 
    {

      g.run([&]{
          this->apply_impl(child0,xData1,yData1,trans);
          this->apply_impl(child1,xData2,yData1,trans);
          });
      g.run_and_wait([&]{
          this->apply_impl(child2,xData1,yData2,trans);
          this->apply_impl(child3,xData2,yData2,trans);
          });
    }
    else
    {
      g.run([&]{
          this->apply_impl(child0,xData1,yData1,trans);
          this->apply_impl(child2,xData2,yData1,trans);
          });
      g.run_and_wait([&]{
          this->apply_impl(child1,xData1,yData2,trans);
          this->apply_impl(child3,xData2,yData2,trans);
          });
    }
  }

//  std::for_each(begin(m_hMatrixData), end(m_hMatrixData),
//                [trans, &X, &Y, this](
//                    const std::pair<shared_ptr<BlockClusterTreeNode<N>>,
//                                    shared_ptr<HMatrixData<ValueType>>>
//                        elem) {
//  
//   IndexRangeType inputRange;
//   IndexRangeType outputRange;
//    if (trans == TransposeMode::NOTRANS) {
//      inputRange = elem.first->data().columnClusterTreeNode->data().indexRange;
//      outputRange = elem.first->data().rowClusterTreeNode->data().indexRange;
//    } else {
//      inputRange = elem.first->data().rowClusterTreeNode->data().indexRange;
//      outputRange = elem.first->data().columnClusterTreeNode->data().indexRange;
//    }
//
//    // Ugly hack since Eigen::Ref constructor does not like const objects
//    Eigen::Ref<Matrix<ValueType>> x_no_const = const_cast<Eigen::Ref<Matrix<ValueType>>&>(X);
//
//    Eigen::Ref<Matrix<ValueType>> xData(x_no_const.block(inputRange[0],0,inputRange[1]-inputRange[0],X.cols()));
//    Eigen::Ref<Matrix<ValueType>> yData(Y.block(outputRange[0],0,outputRange[1]-outputRange[0],Y.cols()));
//  
// 
//    elem.second->apply(xData, yData, trans, 1, 1);
//  
//    //yPermuted.block(outputRange[0],0,outputRange[1]-outputRange[0],yData.cols()) = yData;
//  });
}

template <typename ValueType, int N>
double HMatrix<ValueType, N>::frobeniusNorm_impl(const shared_ptr<BlockClusterTreeNode<N>>& node) const {


   if (node->isLeaf()) return m_hMatrixData.at(node)->frobeniusNorm();


   tbb::task_group g;

   double result = 0;

   double res0;
   double res1;
   double res2;
   double res3;

   g.run([&]{res0 = frobeniusNorm_impl(node->child(0));});
   g.run([&]{res1 = frobeniusNorm_impl(node->child(1));});
   g.run([&]{res2 = frobeniusNorm_impl(node->child(2));});
   g.run_and_wait([&]{res3= frobeniusNorm_impl(node->child(3));});
    
   result = std::sqrt(res0*res0+res1*res1+res2*res2+res3*res3);

   return result;

}


}

#endif
