#ifndef bempp_fmm_octree_hpp
#define bempp_fmm_octree_hpp

#include <complex>
#include <vector>

#include "../common/shared_ptr.hpp"
#include "../fiber/scalar_traits.hpp"
#include "fmm_common.hpp"

#include "../assembly/transposition_mode.hpp"
#include <iostream>
#include <unordered_map>
#include <unordered_set>

namespace Bempp {

class Grid;
}

namespace Fmm {

class Octree {

public:
  // If h is the diameter of the largest element in the grid. Then
  // extend each leaf box on each size by RESIZE_FACTOR * h.
  static constexpr double RESIZE_FACTOR = 1.1;

  // Width of a box relative to the maximum element width h.
  static const int WIDTH_MULTIPLIER = 2.5;

  Octree(const shared_ptr<const Bempp::Grid> &grid, int levels);

  /** \brief Return the number of levels in the grid. */
  unsigned int levels() const;

  /** \brief Return the bounding box of the grid. */
  BoundingBox<double> getBoundingBox() const;

  /** \brief Index of parent box */
  unsigned long getParent(unsigned long nodeIndex) const;

  /** \brief Index of first child */
  unsigned long getFirstChild(unsigned long nodeIndex) const;

  /** \brief Index of last child */
  unsigned long getLastChild(unsigned long nodeIndex) const;

  /** \brief Nodes per side (2^level) */
  unsigned long getNodesPerSide(unsigned int level) const;

  /** \brief Nodes per level (3^level) */
  unsigned long getNodesPerLevel(unsigned int level) const;

  /** \brief Get the number of the octree leaf node that contains the point. */
  unsigned long getLeafContainingPoint(const Point3D<double> &point) const;

  /** \brief Return of a node on a given level is empty. */
  bool isEmpty(unsigned long nodeIndex, unsigned int level) const;

  /** \brief Return the cube width on a given level. */
  double cubeWidth(unsigned int level) const;

  /** \brief Return the extended cube width on a given level. */
  double extendedCubeWidth(unsigned int level) const;

  /** \brief Return the bounds of a specific cube on a given level */
  void cubeBounds(unsigned long nodeIndex, unsigned int level,
                  Vector<double> &lbound, Vector<double> &ubound) const;

  /** \brief Return the extended cube bounds of a specific cube on a given level
   */
  void extendedCubeBounds(unsigned long nodeIndex, unsigned int level,
                          Vector<double> &lbound, Vector<double> &ubound) const;

  const std::vector<unsigned int> &
  getLeafCubeEntities(unsigned long nodeIndex) const;

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

  // Maximum element size in the grid
  double m_maxElementDiam;

  // Store the non empty octree node ids for each level
  std::vector<std::unordered_set<unsigned long>> m_Nodes;

  // Lower bound of global bounding cube
  Vector<double> m_lbound;

  // Upper bounds of global bounding cube
  Vector<double> m_ubound;

  // Stores for each non-empty leaf box the associated grid entities
  std::unordered_map<unsigned long, std::vector<unsigned int>>
      m_leafsToEntities;

  // Value by which the boundaries of boxes are extended to accommodate
  // overhanging triangles. Set to RESIZE_FACTOR * m_maxElementDiam
  double m_extensionSize;
};
}

#include "octree_impl.hpp"

#endif
