// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_BLOCK_CLUSTER_TREE_HPP
#define HMAT_BLOCK_CLUSTER_TREE_HPP

#include "common.hpp"
#include "simple_tree_node.hpp"
#include "bounding_box.hpp"
#include "geometry.hpp"
#include "cluster_tree.hpp"

namespace hmat {

typedef std::function<bool(const BoundingBox &, const BoundingBox &)>
AdmissibilityFunction;

template <int N> struct BlockClusterTreeNodeData {

  BlockClusterTreeNodeData(
      const shared_ptr<const ClusterTreeNode<N>> &rowClusterTreeNode,
      const shared_ptr<const ClusterTreeNode<N>> &columnClusterTreeNode,
      bool admissible);

  shared_ptr<const ClusterTreeNode<N>> rowClusterTreeNode;
  shared_ptr<const ClusterTreeNode<N>> columnClusterTreeNode;

  bool admissible;
};

template <int N >
using BlockClusterTreeNode = SimpleTreeNode<BlockClusterTreeNodeData<N>, N *N>;

typedef BlockClusterTreeNode<2> DefaultBlockClusterTreeNodeType;

template <int N > class BlockClusterTree {

public:
  BlockClusterTree(const shared_ptr<const ClusterTree<N>> &rowClusterTree,
                   const shared_ptr<const ClusterTree<N>> &columnClusterTree,
                   int maxBlockSize,
                   const AdmissibilityFunction &admissibilityFunction);

//  void writeToPdfFile(const std::string &fname, double widthInPoints,
//                      double heightInPoints) const;

  std::size_t rows() const;
  std::size_t columns() const;

  shared_ptr<const BlockClusterTreeNode<N>> root() const;
  shared_ptr<BlockClusterTreeNode<N>> root();

  shared_ptr<const ClusterTree<N>> rowClusterTree() const;
  shared_ptr<const ClusterTree<N>> columnClusterTree() const;

  std::vector<shared_ptr<const BlockClusterTreeNode<N>>> leafNodes() const;
  std::vector<shared_ptr<BlockClusterTreeNode<N>>> leafNodes();

private:
  void initializeBlockClusterTree(
      const AdmissibilityFunction &admissibilityFunction, int maxBlockSize);

  shared_ptr<const ClusterTree<N>> m_rowClusterTree;
  shared_ptr<const ClusterTree<N>> m_columnClusterTree;

  shared_ptr<BlockClusterTreeNode<N>> m_root;
};

template <int N>
void getBlockClusterTreeNodeDimensions(
    const BlockClusterTreeNode<N> &blockClusterTreeNode,
    IndexRangeType &rowClusterRange, IndexRangeType &columnClusterRange,
    std::size_t &numberOfRows, std::size_t &numberOfColumns);

class StandardAdmissibility {
public:
  StandardAdmissibility(double eta);

  bool operator()(const BoundingBox &box1, const BoundingBox &box2) const;

private:
  double m_eta;
};

class WeakAdmissibility {
public:
  bool operator()(const BoundingBox &box1, const BoundingBox &box2) const;
};

typedef BlockClusterTree<2> DefaultBlockClusterTreeType;

}
#include "block_cluster_tree_impl.hpp"

#endif
