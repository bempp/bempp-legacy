#ifndef bempp_fmm_octree_hpp
#define bempp_fmm_octree_hpp

#include <complex>
#include <vector>

#include "../common/shared_ptr.hpp"
#include "../fiber/scalar_traits.hpp"
#include "fmm_common.hpp"

#include "octree_node.hpp"

#include "../assembly/transposition_mode.hpp"
#include <iostream>

namespace Bempp {

class Grid;
}

namespace Fmm {

template <typename CoordinateType> class Octree {

public:
  Octree(const shared_ptr<const Bempp::Grid> &grid, unsigned int levels);

  /** \brief Return the bounding box of the grid. */
  BoundingBox<CoordinateType> getBoundingBox() const;

  /** \brief Index of parent box */
  unsigned long getParent(unsigned long n) const;

  /** \brief Index of first child */
  unsigned long getFirstChild(unsigned long n) const;

  /** \brief Index of last child */
  unsigned long getLastChild(unsigned long n) const;

  /** \brief Nodes per side (2^level) */
  unsigned long getNodesPerSide(unsigned long level) const;

  /** \brief Nodes per level (3^level) */
  unsigned long getNodesPerLevel(unsigned long level) const;

  /** \brief Get the number of the octree leaf node that contains the pont. */
  unsigned long
  getLeafContainingPoint(const Point3D<CoordinateType> &point) const;

private:
  /** \brief return the Morton index of a leaf node */
  unsigned long morton(unsigned long x, unsigned long y, unsigned long z) const;

  /** \brief Overload. */
  unsigned long morton(std::vector<unsigned long> v) const;

  /** \brief Return (x, y, z) indices from Morton index */
  void deMorton(unsigned long *indx, unsigned long *indy, unsigned long *indz,
                unsigned long n) const;

  /** \brief Overload. */
  void deMorton(std::vector<unsigned long> *indv, unsigned long n) const;

  /** \brief Pad an integer with zeros */
  unsigned long dilate3(unsigned long x) const;

  /** \brief Remove padding */
  unsigned long contract3(unsigned long x) const;

  // The underlying grid
  shared_ptr<const Bempp::Grid> m_grid;

  // The number of levels in the Octree
  unsigned int m_levels;

  // Container for Octree nodes
  std::vector<std::vector<OctreeNode<CoordinateType>>> m_OctreeNodes;

  // Global bounding box of grid
  BoundingBox<CoordinateType> m_boundingBox;

  // Size of bounding box across each dimension
  Vector<double> m_boundingBoxSize;

  // Lower bounds of global bounding box
  Vector<double> m_lbound;

  // Upper bounds of global bounding box
  Vector<double> m_ubound;

  // Stores for each leaf box the associated grid entities
  std::vector<std::vector<unsigned int>> m_leafBoxToEntities;

  // Bitfield, which is 1 if a node is empty
  std::vector<std::vector<bool>> m_nodeEmpty;
};
}

#include "octree_impl.hpp"

#endif
