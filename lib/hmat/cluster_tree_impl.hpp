// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_CLUSTER_TREE_IMPL_HPP
#define HMAT_CLUSTER_TREE_IMPL_HPP

#include "common.hpp"
#include "simple_tree_node.hpp"
#include "bounding_box.hpp"
#include "geometry.hpp"
#include "geometry_data_type.hpp"
#include "eigen_fwd.hpp"

#include <CGAL/linear_least_squares_fitting_3.h>

#include "cluster_tree.hpp"

#include <functional>
#include <cassert>

namespace hmat {

inline ClusterTreeNodeData::ClusterTreeNodeData(
    const IndexRangeType &indexRange, const std::vector<Point> &centerPoints, 
      const std::vector<Point>& dofBoundingPoints)
    : indexRange(indexRange) {
   
     geometryData(centerPoints, dofBoundingPoints); 
}

template <int N>
ClusterTree<N>::ClusterTree(const Geometry &geometry, int minBlockSize)
    : m_root(initializeClusterTree(geometry)),
      m_dofPermutation(geometry.size()) {

  splitClusterTreeByGeometry(geometry, m_dofPermutation, minBlockSize);
}

inline void ClusterTreeNodeData::geometryData(const std::vector<Point>& centerPoints,
    const std::vector<Point>& dofBoundingPoints)
{

      Eigen::Matrix<double,3,Eigen::Dynamic> boundingPoints(3,dofBoundingPoints.size());
      
      for (int i = 0; i < dofBoundingPoints.size(); ++i)
      {
        boundingPoints(0,i) = dofBoundingPoints[i].x();
        boundingPoints(1,i) = dofBoundingPoints[i].y();
        boundingPoints(2,i) = dofBoundingPoints[i].z();
      }
  
      linear_least_squares_fitting_3(centerPoints.begin(),
        centerPoints.end(),mainLine,centroid,CGAL::Dimension_tag<0>());

      // Now compute the bounding box
      
      double xmin, xmax, ymin, ymax, zmin, zmax;

      Eigen::Vector3d maxVector = boundingPoints.rowwise().maxCoeff();
      Eigen::Vector3d minVector = boundingPoints.rowwise().minCoeff();

      xmin = minVector(0); ymin = minVector(1); zmin = minVector(2);
      xmax = maxVector(0); ymax = maxVector(1); zmax = maxVector(2);

      boundingBox = BoundingBox(xmin,xmax,ymin,ymax,zmin,zmax);

      // Store the diameter
      
      diameter = clusterDiameter(dofBoundingPoints);

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

  std::vector<Point> centerPoints;
  centerPoints.reserve(geometry.size());
  std::vector<Point> dofBoundingPoints;
  dofBoundingPoints.reserve(8*geometry.size());

  for (const auto &geometryData : geometry)
  {
    centerPoints.push_back(geometryData->center);
    geometryData->boundingBox.corners(dofBoundingPoints);
  }

  return make_shared<ClusterTreeNode<N>>(ClusterTreeNodeData(indexRange,
       centerPoints,dofBoundingPoints));
}

template <int N>
const std::vector<std::size_t> &
ClusterTree<N>::hMatDofToOriginalDofMap() const {
  return m_dofPermutation.hMatDofToOriginalDofMap();
}

template <int N>
const std::vector<std::size_t> &
ClusterTree<N>::originalDofToHMatDofMap() const {
  return m_dofPermutation.originalDofToHMatDofMap();
}

template <>
inline void
ClusterTree<2>::splitClusterTreeByGeometry(const Geometry &geometry,
                                           DofPermutation &dofPermutation,
                                           int minBlockSize) {

  std::function<void(const shared_ptr<ClusterTreeNode<2>> & clusterTreeNode,
                     const IndexSetType & indexSet)> splittingFun;

  splittingFun = [&dofPermutation, &geometry, minBlockSize, &splittingFun](
      const shared_ptr<ClusterTreeNode<2>> &clusterTreeNode,
      const IndexSetType &indexSet) {

    std::size_t indexSetSize = indexSet.size();

    assert(indexSetSize == clusterTreeNode->data().indexRange[1] -
                               clusterTreeNode->data().indexRange[0]);

    if (indexSetSize > minBlockSize) {

      IndexSetType firstIndexSet;
      IndexSetType secondIndexSet;

      // Approximate the points in the cluster with a line.
      
      Line line;
      Point centroid;

      auto plane = clusterTreeNode->data().mainLine.perpendicular_plane(
          clusterTreeNode->data().centroid);
      std::vector<Point> firstPointSet;
      std::vector<Point> secondPointSet;
      std::vector<Point> firstBoundingPointSet;
      std::vector<Point> secondBoundingPointSet;

      for (auto index: indexSet)
       if (plane.oriented_side(geometry[index]->center)==CGAL::ON_POSITIVE_SIDE)
       {
         firstIndexSet.push_back(index);
         firstPointSet.push_back(geometry[index]->center);
         geometry[index]->boundingBox.corners(firstBoundingPointSet);
       }
     else
       {
         secondIndexSet.push_back(index);
         secondPointSet.push_back(geometry[index]->center);
         geometry[index]->boundingBox.corners(secondBoundingPointSet); 
       }

      if (firstIndexSet.size()==0) throw std::runtime_error(
          "hmat::splitClusterTreeByGeometry(): First index set is zero.");

      if (secondIndexSet.size()==0) throw std::runtime_error(
          "hmat::splitClusterTreeByGeometry(): Second index set is zero.");

      auto pivot = firstIndexSet.size();

      IndexRangeType newRangeFirst = clusterTreeNode->data().indexRange;
      IndexRangeType newRangeSecond = clusterTreeNode->data().indexRange;

      newRangeFirst[1] = newRangeSecond[0] = newRangeFirst[0] + pivot;

      clusterTreeNode->addChild(
          ClusterTreeNodeData(newRangeFirst, firstPointSet,
            firstBoundingPointSet), 0);
      clusterTreeNode->addChild(
          ClusterTreeNodeData(newRangeSecond, secondPointSet,
            secondBoundingPointSet), 1);
      splittingFun(clusterTreeNode->child(0), firstIndexSet);
      splittingFun(clusterTreeNode->child(1), secondIndexSet);

    } else {

      int originalIndexCount = 0;
      const auto &indexRange = clusterTreeNode->data().indexRange;
      for (int hMatDof = indexRange[0]; hMatDof < indexRange[1]; ++hMatDof) {
        dofPermutation.addDofIndexPair(indexSet[originalIndexCount], hMatDof);
        ++originalIndexCount;
      }
    }

  };
  splittingFun(m_root, fillIndexRange(0, geometry.size()));
  std::cout << "Cluster tree generated." << std::endl;
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
