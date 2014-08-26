// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_CLUSTER_TREE_IMPL_HPP
#define HMAT_CLUSTER_TREE_IMPL_HPP

#include "common.hpp"
#include "simple_tree_node.hpp"
#include "bounding_box.hpp"
#include "geometry.hpp"
#include "geometry_data_type.hpp"

#include "cluster_tree.hpp"

#include <functional>

namespace hmat {

inline ClusterTreeNodeData::ClusterTreeNodeData(
    const IndexRangeType &indexRange, const BoundingBox &boundingBox)
    : indexRange(indexRange), boundingBox(boundingBox) {}

template <int N>
ClusterTree<N>::ClusterTree(const Geometry &geometry, int maxSize)
    : m_root(initializeClusterTree(geometry)),
      m_dofPermutation(geometry.size()) {

  splitClusterTreeByGeometry(geometry, m_dofPermutation, maxSize);
}

template <int N> std::size_t ClusterTree<N>::numberOfDofs() const {
  return (m_root->data().indexRange[1] - m_root->data().indexRange[0]);
}

template <int N>
const shared_ptr<const ClusterTreeNode<N>> ClusterTree<N>::root() const {
  return m_root;
}

template <int N> const shared_ptr<ClusterTreeNode<N>> ClusterTree<N>::root() {
  return m_root;
}

template <int N>
shared_ptr<ClusterTreeNode<N>>
ClusterTree<N>::initializeClusterTree(const Geometry &geometry) {

  IndexRangeType indexRange{{0, geometry.size()}};

  BoundingBox b;

  for (const auto &geometryData : geometry)
    b.merge(geometryData->boundingBox);

  return make_shared<ClusterTreeNode<N>>(ClusterTreeNodeData(indexRange, b));
}

template <int N>
const std::vector<std::size_t>& ClusterTree<N>::hMatDofToOriginalDofMap() const {
  return m_dofPermutation.hMatDofToOriginalDofMap();
}

template <int N>
const std::vector<std::size_t>& ClusterTree<N>::originalDofToHMatDofMap() const {
  return m_dofPermutation.originalDofToHMatDofMap();
}


template <>
void ClusterTree<2>::splitClusterTreeByGeometry(const Geometry &geometry,
                                                DofPermutation &dofPermutation,
                                                int maxSize) {

  std::function<void(const shared_ptr<ClusterTreeNode<2>> & clusterTreeNode,
                     const IndexSetType & indexSet)> splittingFun;

  splittingFun = [&dofPermutation, &geometry, maxSize, &splittingFun](
      const shared_ptr<ClusterTreeNode<2>> &clusterTreeNode,
      const IndexSetType &indexSet) {

    if (indexSet.size() > maxSize) {

      IndexSetType firstIndexSet;
      IndexSetType secondIndexSet;

      auto boxes = clusterTreeNode->data().boundingBox.divideMaxDimension(.5);

      auto firstBoundingBox = boxes.first;
      auto secondBoundingBox = boxes.second;

      for (auto index : indexSet)
        if (firstBoundingBox.contains(geometry[index]->center))
          firstIndexSet.push_back(index);
        else
          secondIndexSet.push_back(index);

      if (secondIndexSet.size() == 0) {
        clusterTreeNode->data().boundingBox = firstBoundingBox;
        splittingFun(clusterTreeNode, indexSet);
        return;
      }
      if (firstIndexSet.size() == 0) {
        clusterTreeNode->data().boundingBox = secondBoundingBox;
        splittingFun(clusterTreeNode, indexSet);
        return;
      }

      auto pivot = firstIndexSet.size();

      IndexRangeType newRangeFirst = clusterTreeNode->data().indexRange;
      IndexRangeType newRangeSecond = clusterTreeNode->data().indexRange;

      newRangeFirst[1] = newRangeFirst[0] + pivot;
      newRangeSecond[0] = newRangeFirst[0] + pivot;

      clusterTreeNode->addChild(
          ClusterTreeNodeData(newRangeFirst, firstBoundingBox), 0);
      clusterTreeNode->addChild(
          ClusterTreeNodeData(newRangeSecond, secondBoundingBox), 1);
      splittingFun(clusterTreeNode->child(0), firstIndexSet);
      splittingFun(clusterTreeNode->child(1), secondIndexSet);

      clusterTreeNode->data().boundingBox =
          clusterTreeNode->child(0)->data().boundingBox;
      clusterTreeNode->data().boundingBox.merge(
          clusterTreeNode->child(1)->data().boundingBox);
    } else {

      BoundingBox b;
      for (const auto &elem : indexSet)
        b.merge(geometry[elem]->boundingBox);

      int originalIndexCount = 0;
      const auto &indexRange = clusterTreeNode->data().indexRange;
      for (int hMatDof = indexRange[0]; hMatDof < indexRange[1]; ++hMatDof) {
        dofPermutation.addDofIndexPair(indexSet[originalIndexCount], hMatDof);
        ++originalIndexCount;
      }
    }
  };
  splittingFun(m_root, fillIndexRange(0, geometry.size()));
}

template <int N>
std::size_t
ClusterTree<N>::mapOriginalDofToHMatDof(std::size_t originalDofIndex) const {
  return m_dofPermutation.mapOriginalDofToHMatDof(originalDofIndex);
}

template <int N>
std::size_t
ClusterTree<N>::mapHMatDofToOriginalDof(std::size_t hMatDofIndex) const {
  return m_dofPermutation.mapHMatDofToOriginalDof(hMatDofIndex);
}

template <int N>
std::vector<shared_ptr<const ClusterTreeNode<N>>>
ClusterTree<N>::leafNodes() const {

  return m_root->leafNodes;
}

template <int N>
std::vector<shared_ptr<ClusterTreeNode<N>>> ClusterTree<N>::leafNodes() {

  return m_root->leafNodes();
}
}
#endif
