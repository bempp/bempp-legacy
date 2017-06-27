#ifndef bempp_fmm_octree_impl_hpp
#define bempp_fmm_octree_impl_hpp

#include "../grid/grid.hpp"
#include "./octree.hpp"

namespace Fmm {

inline Octree::Octree(const shared_ptr<const Bempp::Grid> &grid,
                      unsigned int levels)
    : m_grid(grid), m_levels(levels) {

  // Compute maximum allowed level

  auto gridView = m_grid->leafView();

  // Compute the grid bounding box.
  Vector<double> lbound(3);
  Vector<double> ubound(3);

  m_grid->getBoundingBox(m_lbound, m_ubound);

  // Compute width

  double width = (m_ubound - m_lbound).maxCoeff();

  // Extend ubound to the maximum width in each dimension

  m_ubound = m_lbound + width * Vector<double>::Ones(3);

  // Get the maximum allowed level

  m_maxElementDiam = gridView->maximumElementDiameter();
  m_extensionSize = m_maxElementDiam * Octree::RESIZE_FACTOR;

  int maxLevels = (int)std::trunc(
      std::log2(width / (Octree::WIDTH_MULTIPLIER * m_maxElementDiam)));

  if (levels < 0 || levels > maxLevels)
    m_levels = maxLevels;

  // Compute node assignments

  m_Nodes.resize(m_levels);

  auto view = grid->leafView();
  const auto &indexSet = view->indexSet();

  for (auto it = view->entityIterator<0>(); !it->finished(); it->next()) {
    const auto &entity = it->entity();
    int entityIndex = indexSet.entityIndex(entity);

    Vector<double> centroid(3);

    // Get the center and corners of the entity.
    entity.geometry().getCenter(centroid);

    // Find out in which leaf node the element lives
    unsigned long nodeIndex =
        getLeafContainingPoint({centroid(0), centroid(1), centroid(2)});

    // Check if box is empty.
    if (m_Nodes[m_levels - 1].count(nodeIndex) == 0) {
      // Box has not yet been inserted as nonzero
      m_Nodes[m_levels - 1].insert(nodeIndex);
      // Mark all parents as nonzero
      unsigned long parent = nodeIndex;
      for (int level = m_levels; level >= 1; level--) {
        parent = getParent(parent);
        m_Nodes[level - 1].insert(parent);
      }
    }

    // Add entity index to the list of entities associated with the leaf box.
    if (m_leafsToEntities.count(nodeIndex) == 0)
      m_leafsToEntities[nodeIndex] == std::vector<unsigned int>();

    m_leafsToEntities[nodeIndex].push_back(entityIndex);
  }
}

inline BoundingBox<double> Octree::getBoundingBox() const {

  return Bempp::createBoundingBox(m_lbound(0), m_lbound(1), m_lbound(2),
                                  m_ubound(0), m_ubound(1), m_ubound(2));
}

inline unsigned int Octree::levels() const { return m_levels; }

inline unsigned long Octree::getParent(unsigned long nodeIndex) const {
  return nodeIndex >> 3;
}

inline unsigned long Octree::getFirstChild(unsigned long nodeIndex) const {
  return nodeIndex << 3;
}

inline unsigned long Octree::getLastChild(unsigned long nodeIndex) const {
  return (nodeIndex << 3) + 7;
}

inline unsigned long Octree::getNodesPerSide(unsigned int level) const {
  return 1 << level; // 2^level
}

inline unsigned long Octree::getNodesPerLevel(unsigned int level) const {
  return 1 << 3 * level; // 8^level;
}

inline unsigned long
Octree::getLeafContainingPoint(const Point3D<double> &point) const {
  int invleafsize = getNodesPerSide(m_levels);

  Vector<double> boxSize = m_ubound - m_lbound;

  // be careful of precision, outside allocation bad
  Vector<double> pt(3);
  pt(0) = (point.x - m_lbound(0)) / boxSize(0);
  pt(1) = (point.y - m_lbound(1)) / boxSize(1);
  pt(2) = (point.z - m_lbound(2)) / boxSize(2);

  double zero = 0;

  std::vector<unsigned long> ind;
  for (int i = 0; i < 3; ++i)
    ind.push_back(
        std::min(int(std::max(zero, pt(i)) * invleafsize), invleafsize - 1));

  return morton(ind);
}

inline bool Octree::isEmpty(unsigned long nodeIndex, unsigned int level) const {

  return (m_Nodes[level - 1].count(nodeIndex) == 0);
}

inline double Octree::cubeWidth(unsigned int level) const {

  double width = m_ubound(0) - m_lbound(0);

  return width / (1 << level);
}

inline double Octree::extendedCubeWidth(unsigned int level) const {

  double width = m_ubound(0) - m_lbound(0);

  return width / (1 << level) + 2 * m_extensionSize;
}

inline void Octree::cubeBounds(unsigned long nodeIndex, unsigned int level,
                               Vector<double> &lbound,
                               Vector<double> &ubound) const {

  unsigned long indx, indy, indz;
  deMorton(&indx, &indy, &indz, nodeIndex);

  double width = cubeWidth(level);

  lbound.resize(3);
  ubound.resize(3);
  lbound(0) = m_lbound(0) + indx * width;
  lbound(1) = m_lbound(1) + indy * width;
  lbound(2) = m_lbound(2) + indz * width;

  ubound(0) = m_lbound(0) + (indx + 1) * width;
  ubound(1) = m_lbound(1) + (indy + 1) * width;
  ubound(2) = m_lbound(2) + (indz + 1) * width;
}

inline void Octree::extendedCubeBounds(unsigned long nodeIndex,
                                       unsigned int level,
                                       Vector<double> &lbound,
                                       Vector<double> &ubound) const {

  cubeBounds(nodeIndex, level, lbound, ubound);

  lbound.resize(3);
  ubound.resize(3);

  lbound(0) -= m_extensionSize;
  lbound(1) -= m_extensionSize;
  lbound(2) -= m_extensionSize;

  ubound(0) += m_extensionSize;
  ubound(1) += m_extensionSize;
  ubound(2) += m_extensionSize;
}

// template <typename CoordinateType> double cubeWidth(unsigned int level)
// const
// {}

// Dilate an integer, in between each and every bit of the number inserting
// two zero bits
inline unsigned long Octree::dilate3(unsigned long x) const {
  if (x > 0x000003FF)
    throw std::invalid_argument("dilate3(x): argument x"
                                "exceeds maximum allowed (1023)");
  // x = ---- ---- ---- ---- ---- --98 7654 3210
  x = (x | (x << 16)) &
      0x030000FF; // x = ---- --98 ---- ---- ---- ---- 7654 3210
  x = (x | (x << 8)) &
      0x0300F00F; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
  x = (x | (x << 4)) &
      0x030C30C3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
  x = (x | (x << 2)) &
      0x09249249; // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0

  return x;
}

// undo Dilate, trashing padding bits
inline unsigned long Octree::contract3(unsigned long x) const {
  x &= 0x09249249; // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
  x = (x | (x >> 2)) &
      0x030C30C3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
  x = (x | (x >> 4)) &
      0x0300F00F; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
  x = (x | (x >> 8)) &
      0x030000FF; // x = ---- --98 ---- ---- ---- ---- 7654 3210
  x = (x | (x >> 16)) &
      0x000003FF; // x = ---- ---- ---- ---- ---- --98 7654 3210

  return x;
}

inline unsigned long Octree::morton(unsigned long x, unsigned long y,
                                    unsigned long z) const {
  return dilate3(x) | (dilate3(y) << 1) | (dilate3(z) << 2);
}

inline unsigned long Octree::morton(std::vector<unsigned long> v) const {
  return morton(v[0], v[1], v[2]);
}

inline void Octree::deMorton(unsigned long *indx, unsigned long *indy,
                             unsigned long *indz, unsigned long n) const {
  *indx = contract3(n);
  *indy = contract3(n >> 1);
  *indz = contract3(n >> 2);
}
}

#endif
