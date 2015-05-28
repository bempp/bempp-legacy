// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_CLUSTER_TREE_HPP
#define HMAT_CLUSTER_TREE_HPP

#include "common.hpp"
#include "simple_tree_node.hpp"
#include "geometry.hpp"
#include "dof_permutation.hpp"

namespace hmat {

struct ClusterTreeNodeData {

  ClusterTreeNodeData(const IndexRangeType &indexRange,
                      const std::vector<Point>& clusterPoints);

  IndexRangeType indexRange;
  std::vector<Point> clusterPoints;
  double diameter;
};

template <int N> using ClusterTreeNode = SimpleTreeNode<ClusterTreeNodeData, N>;

typedef ClusterTreeNode<2> DefaultClusterTreeNodeType;

template <int N> class ClusterTree {

public:
  ClusterTree(const Geometry &geometry, int minBlockSize);

  const shared_ptr<const ClusterTreeNode<N>> root() const;
  const shared_ptr<ClusterTreeNode<N>> root();

  std::size_t mapOriginalDofToHMatDof(std::size_t originalDofIndex) const;
  std::size_t mapHMatDofToOriginalDof(std::size_t hMatDofIndex) const;

  const std::vector<std::size_t> &hMatDofToOriginalDofMap() const;
  const std::vector<std::size_t> &originalDofToHMatDofMap() const;

  std::vector<shared_ptr<const ClusterTreeNode<N>>> leafNodes() const;
  std::vector<shared_ptr<ClusterTreeNode<N>>> leafNodes();

  std::size_t numberOfDofs() const;

private:
  shared_ptr<ClusterTreeNode<N>>
  initializeClusterTree(const Geometry &geometry);
  void splitClusterTreeByGeometry(const Geometry &geometry,
                                  DofPermutation &dofPermutation,
                                  int minBlockSize);

  shared_ptr<ClusterTreeNode<N>> m_root;
  DofPermutation m_dofPermutation;
};

typedef ClusterTree<2> DefaultClusterTreeType;
}
#include "cluster_tree_impl.hpp"

#endif
