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
