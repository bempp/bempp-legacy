#ifndef bempp_fmm_octree_impl_hpp
#define bempp_fmm_octree_impl_hpp

#include "./octree.hpp"

namespace Fmm {

template <typename CoordinateType>
Octree<CoordinateType>::Octree(const shared_ptr<const Bempp::Grid>& grid, unsigned int levels)
    : m_grid(grid)
    , m_levels(levels)
{
    m_grid->getBoundingBox(m_boundingBox);
    m_boundingBoxSize = Bempp::getBoundingBoxSize(m_boundingBox);
    m_lbound.resize(3);
    m_ubound.resize(3);

    m_lbound(0) = m_boundingBox.lbound.x;
    m_lbound(1) = m_boundingBox.lbound.y;
    m_lbound(2) = m_boundingBox.lbound.z;

    m_ubound(0) = m_boundingBox.ubound.x;
    m_ubound(1) = m_boundingBox.ubound.y;
    m_ubound(2) = m_boundingBox.ubound.z;

    m_OctreeNodes.resize(m_levels);

    auto boxSize = getBoundingBoxSize(m_boundingBox);

    // initialise octree stucture (don't bother storing the lowest two levels)
    for (unsigned int level = 1; level <= m_levels; ++level) {
        unsigned int nNodesOnSide = getNodesPerSide(level);
        unsigned int nNodes = getNodesPerLevel(level);
        m_OctreeNodes[level - 1].reserve(nNodes);
        for (unsigned int nodeNumber = 0; nodeNumber < nNodes; ++nodeNumber) {
            unsigned long ind[3];
            deMorton(&ind[0], &ind[1], &ind[2], nodeNumber);

            double xmin = ind[0] * boxSize(0) + m_boundingBox.lbound.x;
            double ymin = ind[1] * boxSize(1) + m_boundingBox.lbound.y;
            double zmin = ind[2] * boxSize(2) + m_boundingBox.lbound.z;

            double xmax = (ind[0] + 1) * boxSize(0) + m_boundingBox.lbound.x;
            double ymax = (ind[1] + 1) * boxSize(1) + m_boundingBox.lbound.y;
            double zmax = (ind[2] + 1) * boxSize(2) + m_boundingBox.lbound.z;

            m_OctreeNodes[level - 1].push_back(
                OctreeNode<CoordinateType>(nodeNumber, level,
                    Bempp::createBoundingBox<CoordinateType>(xmin, ymin, zmin, xmax, ymax, zmax)));
        }
    }
}

template <typename CoordinateType>
BoundingBox<CoordinateType> Octree<CoordinateType>::getBoundingBox() const
{

    return m_boundingBox;
}

template <typename CoordinateType>
unsigned long Octree<CoordinateType>::getParent(unsigned long n) const
{
    return n >> 3;
}

template <typename CoordinateType>
unsigned long Octree<CoordinateType>::getFirstChild(unsigned long n) const
{
    return n << 3;
}

template <typename CoordinateType>
unsigned long Octree<CoordinateType>::getLastChild(unsigned long n) const
{
    return (n << 3) + 7;
}

template <typename CoordinateType>
unsigned long Octree<CoordinateType>::getNodesPerSide(unsigned long level) const
{
    return 1 << level; // 2^level
}

template <typename CoordinateType>
unsigned long Octree<CoordinateType>::getNodesPerLevel(unsigned long level) const
{
    return 1 << 3 * level; // 8^level;
}

template <typename CoordinateType>
unsigned long Octree<CoordinateType>::getLeafContainingPoint(
    const Point3D<CoordinateType>& point) const
{
    int invleafsize = getNodesPerSide(m_levels);

    // be careful of precision, outside allocation bad
    Vector<CoordinateType> pt(3);
    for (int i = 0; i < 3; ++i)
        pt(i) = (point(i) - m_lbound(i)) / m_boundingBoxSize(i);

    CoordinateType zero = CoordinateType(0);

    std::vector<unsigned long> ind;
    for (int i = 0; i < 3; ++i)
        ind.push_back(std::min(int(std::max(zero, pt(i)) * invleafsize),
            invleafsize - 1));

    return morton(ind);
}

// Dilate an integer, in between each and every bit of the number inserting
// two zero bits
template <typename CoordinateType>
unsigned long Octree<CoordinateType>::dilate3(unsigned long x) const
{
    if (x > 0x000003FF)
        throw std::invalid_argument("dilate3(x): argument x"
                                    "exceeds maximum allowed (1023)");
    // x = ---- ---- ---- ---- ---- --98 7654 3210
    x = (x | (x << 16)) & 0x030000FF; // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x | (x << 8)) & 0x0300F00F;  // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x | (x << 4)) & 0x030C30C3;  // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x | (x << 2)) & 0x09249249;  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0

    return x;
}

// undo Dilate, trashing padding bits
template <typename CoordinateType>
unsigned long Octree<CoordinateType>::contract3(unsigned long x) const
{
    x &= 0x09249249;                  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    x = (x | (x >> 2)) & 0x030C30C3;  // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x | (x >> 4)) & 0x0300F00F;  // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x | (x >> 8)) & 0x030000FF;  // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x | (x >> 16)) & 0x000003FF; // x = ---- ---- ---- ---- ---- --98 7654 3210

    return x;
}

template <typename CoordinateType>
unsigned long Octree<CoordinateType>::morton(unsigned long x, unsigned long y, unsigned long z) const
{
    return dilate3(x) | (dilate3(y) << 1) | (dilate3(z) << 2);
}

template <typename CoordinateType>
unsigned long Octree<CoordinateType>::morton(std::vector<unsigned long> v) const
{
    return morton(v[0], v[1], v[2]);
}

template <typename CoordinateType>
void Octree<CoordinateType>::deMorton(unsigned long* indx, unsigned long* indy,
    unsigned long* indz, unsigned long n) const
{
    *indx = contract3(n);
    *indy = contract3(n >> 1);
    *indz = contract3(n >> 2);
}
}

#endif
