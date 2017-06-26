#ifndef bempp_fmm_octree_node_hpp
#define bempp_fmm_octree_node_hpp

#include <complex>
#include <vector>

#include "fmm_common.hpp"

#include "../common/shared_ptr.hpp"
#include "../fiber/scalar_traits.hpp"

namespace Fmm {

/** \cond FORWARD_DECL */
template <typename ResultType> class Octree;
/** \endcond */

template <typename CoordinateType> class OctreeNode {

  friend class Octree<CoordinateType>;

public:
  OctreeNode(unsigned long number, unsigned int level,
             const BoundingBox<CoordinateType> &boundingBox,
             const Octree<CoordinateType> &octree);

private:
  unsigned long m_number;
  unsigned long m_level;
  BoundingBox<CoordinateType> m_boundingBox;
  const Octree<CoordinateType> &m_octree;
};
}
#include "./octree_node_impl.hpp"

#endif
