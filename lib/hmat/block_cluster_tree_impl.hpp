// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_BLOCK_CLUSTER_TREE_IMPL_HPP
#define HMAT_BLOCK_CLUSTER_TREE_IMPL_HPP

#include "block_cluster_tree.hpp"
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

//#include "cairo/cairo.h"
//#include "cairo/cairo-pdf.h"

namespace hmat {

template <int N>
BlockClusterTreeNodeData<N>::BlockClusterTreeNodeData(
    const shared_ptr<const ClusterTreeNode<N>> &rowClusterTreeNode,
    const shared_ptr<const ClusterTreeNode<N>> &columnClusterTreeNode,
    bool admissible)
    : rowClusterTreeNode(rowClusterTreeNode),
      columnClusterTreeNode(columnClusterTreeNode), admissible(admissible) {}

template <int N>
BlockClusterTree<N>::BlockClusterTree(
    const shared_ptr<const ClusterTree<N>> &rowClusterTree,
    const shared_ptr<const ClusterTree<N>> &columnClusterTree, int maxBlockSize,
    const AdmissibilityFunction &admissibilityFunction)
    : m_rowClusterTree(rowClusterTree), m_columnClusterTree(columnClusterTree) {

  initializeBlockClusterTree(admissibilityFunction, maxBlockSize);
}

// template <int N>
// void BlockClusterTree<N>::writeToPdfFile(const std::string &fname,
//                                         double widthInPoints,
//                                         double heightInPoints) const {
//
//  cairo_surface_t *surface;
//  surface =
//      cairo_pdf_surface_create(fname.c_str(), widthInPoints, heightInPoints);
//  cairo_t *cr = cairo_create(surface);
//  cairo_set_line_width(cr, 1);
//
//  std::function<void(shared_ptr<BlockClusterTreeNode<N>>)> paintImpl;
//
//  paintImpl = [&cr, &paintImpl, widthInPoints, heightInPoints, this](
//      shared_ptr<BlockClusterTreeNode<N>> node) {
//
//    auto nodeData = node->data();
//
//    if (node->isLeaf()) {
//
//      auto rowIndexRange = nodeData.rowClusterTreeNode->data().indexRange;
//      auto columnIndexRange =
//      nodeData.columnClusterTreeNode->data().indexRange;
//
//      auto rowDiameter = rowIndexRange[1] - rowIndexRange[0];
//      auto columnDiameter = columnIndexRange[1] - columnIndexRange[0];
//
//      double x = (widthInPoints * columnIndexRange[0]) / (this->columns());
//      double y = (heightInPoints * rowIndexRange[0]) / (this->rows());
//
//      double rectWidth = (widthInPoints * columnDiameter) / this->columns();
//      double rectHeight = (heightInPoints * rowDiameter) / this->rows();
//
//      double green = nodeData.admissible ? 1.0 : 0.0;
//      double red = 1 - green;
//
//      cairo_set_source_rgb(cr, red, green, 0);
//      cairo_rectangle(cr, x, y, rectWidth, rectHeight);
//      cairo_fill_preserve(cr);
//      cairo_set_source_rgb(cr, 0, 0, 0);
//      cairo_stroke(cr);
//    } else {
//      for (int i = 0; i < N * N; ++i)
//        paintImpl(node->child(i));
//    }
//  };
//
//  paintImpl(m_root);
//
//  cairo_destroy(cr);
//  cairo_surface_destroy(surface);
//}

template <int N> std::size_t BlockClusterTree<N>::rows() const {
  return m_rowClusterTree->numberOfDofs();
}

template <int N> std::size_t BlockClusterTree<N>::columns() const {
  return m_columnClusterTree->numberOfDofs();
}

template <int N>
shared_ptr<const BlockClusterTreeNode<N>> BlockClusterTree<N>::root() const {
  return m_root;
}

template <int N>
shared_ptr<BlockClusterTreeNode<N>> BlockClusterTree<N>::root() {
  return m_root;
}

template <int N>
shared_ptr<const ClusterTree<N>> BlockClusterTree<N>::rowClusterTree() const {

  return m_rowClusterTree;
}

template <int N>
shared_ptr<const ClusterTree<N>>
BlockClusterTree<N>::columnClusterTree() const {
  return m_columnClusterTree;
}

template <int N>
std::vector<shared_ptr<const BlockClusterTreeNode<N>>>
BlockClusterTree<N>::leafNodes() const {
  return m_root->leafNodes();
}

template <int N>
std::vector<shared_ptr<BlockClusterTreeNode<N>>>
BlockClusterTree<N>::leafNodes() {
  return m_root->leafNodes();
}

template <int N>
void BlockClusterTree<N>::initializeBlockClusterTree(
    const AdmissibilityFunction &admissibilityFunction, int maxBlockSize) {

  std::function<void(const shared_ptr<BlockClusterTreeNode<N>> &)>
      splittingFunction;

  splittingFunction =
      [this, &admissibilityFunction, &splittingFunction, &maxBlockSize](
          const shared_ptr<BlockClusterTreeNode<N>> &node) {

        auto nodeData = node->data();

        // Adjust admissibility condition to only accept blocks smaller than
        // maxBlockSize

        auto rowClusterTreeNodeIndexRange =
            nodeData.rowClusterTreeNode->data().indexRange;
        auto columnClusterTreeNodeIndexRange =
            nodeData.columnClusterTreeNode->data().indexRange;
        auto rowBlockSize =
            rowClusterTreeNodeIndexRange[1] - rowClusterTreeNodeIndexRange[0];
        auto columnBlockSize = columnClusterTreeNodeIndexRange[1] -
                               columnClusterTreeNodeIndexRange[0];

        if (columnBlockSize > maxBlockSize || rowBlockSize > maxBlockSize)
          nodeData.admissible = false;

        // If admissible do not refine further

        if (nodeData.admissible)
          return;

        auto rowClusterData = nodeData.rowClusterTreeNode->data();
        auto columnClusterData = nodeData.columnClusterTreeNode->data();

        // If row or column cluster is leaf do not refine further

        if (nodeData.rowClusterTreeNode->isLeaf() ||
            nodeData.columnClusterTreeNode->isLeaf())
          return;

        // Create the block clusters

        for (int rowCount = 0; rowCount < N; ++rowCount) {
          auto rowChild = nodeData.rowClusterTreeNode->child(rowCount);
          for (int columnCount = 0; columnCount < N; ++columnCount) {
            auto columnChild =
                nodeData.columnClusterTreeNode->child(columnCount);
            node->addChild(
                BlockClusterTreeNodeData<N>(
                    rowChild, columnChild,
                    admissibilityFunction(rowClusterData, columnClusterData)),
                N * rowCount + columnCount);
            splittingFunction(node->child(N * rowCount + columnCount));
          }
        }
        // tbb::parallel_for(tbb::blocked_range<int>(0,N*N),[&splittingFunction,&node](
        //      const tbb::blocked_range<int>& range){
        //    for (int i = range.begin(); i!=range.end(); ++i)
        //      splittingFunction(node->child(i));});
      };

  bool admissible = admissibilityFunction(m_rowClusterTree->root()->data(),
                                          m_columnClusterTree->root()->data());
  m_root = shared_ptr<BlockClusterTreeNode<N>>(
      new BlockClusterTreeNode<N>(BlockClusterTreeNodeData<N>(
          m_rowClusterTree->root(), m_columnClusterTree->root(), admissible)));
  splittingFunction(m_root);
  std::cout << "BlockClusterTree generated." << std::endl;
}

template <int N>
void getBlockClusterTreeNodeDimensions(
    const BlockClusterTreeNode<N> &blockClusterTreeNode,
    IndexRangeType &rowClusterRange, IndexRangeType &columnClusterRange,
    std::size_t &numberOfRows, std::size_t &numberOfColumns) {

  rowClusterRange =
      blockClusterTreeNode.data().rowClusterTreeNode->data().indexRange;
  columnClusterRange =
      blockClusterTreeNode.data().columnClusterTreeNode->data().indexRange;
  numberOfRows = rowClusterRange[1] - rowClusterRange[0];
  numberOfColumns = columnClusterRange[1] - columnClusterRange[0];
}

inline StrongAdmissibility::StrongAdmissibility(double eta) : m_eta(eta) {}

inline bool StrongAdmissibility::
operator()(const ClusterTreeNodeData &cluster1,
           const ClusterTreeNodeData &cluster2) const {
  double diam1 = cluster1.diameter;
  double diam2 = cluster2.diameter;

  double dist = cluster1.boundingBox.distance(cluster2.boundingBox);

  return std::min(diam1, diam2) < m_eta * dist;
}

inline bool WeakAdmissibility::
operator()(const ClusterTreeNodeData &cluster1,
           const ClusterTreeNodeData &cluster2) const {

  return cluster1.boundingBox.distance(cluster2.boundingBox) > 0;
}

inline HighFrequencyAdmissibility::HighFrequencyAdmissibility(double k): m_k(k){}

inline bool HighFrequencyAdmissibility::
operator()(const ClusterTreeNodeData &cluster1,
           const ClusterTreeNodeData &cluster2) const {

  double dist = cluster1.boundingBox.distance(cluster2.boundingBox);
  double diam1 = cluster1.diameter;
  double diam2 = cluster2.diameter;
  return m_k*diam1*diam2/dist < 10; 
}

}
#endif
