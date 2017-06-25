#ifndef bempp_fmm_octree_node_impl_hpp
#define bempp_fmm_octree_node_impl_hpp

#include "./octree_node.hpp"

namespace Fmm {

template <typename CoordinateType>
OctreeNode<CoordinateType>::OctreeNode(unsigned long number, unsigned int level)
    : m_number(number)
    , m_level(level)
{
}
}

#endif
