#ifndef bempp_fmm_octree_hpp
#define bempp_fmm_octree_hpp

#include <complex>
#include <vector>

#include "../common/shared_ptr.hpp"
#include "../fiber/scalar_traits.hpp"
#include "fmm_common.hpp"

#include "dof_permutation.hpp"
#include "octree_node.hpp"

#include "../assembly/transposition_mode.hpp"
#include <iostream>

namespace Bempp {

class Grid;
}

namespace Fmm {

template <typename CoordinateType>
class Octree {

public:
    Octree(const shared_ptr<const Bempp::Grid>& grid, unsigned int levels);

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

private:
    /** \brief return the Morton index of a leaf node */
    unsigned long morton(unsigned long x, unsigned long y, unsigned long z) const;

    /** \brief Overload. */
    unsigned long morton(std::vector<unsigned long> v) const;

    /** \brief Return (x, y, z) indices from Morton index */
    void deMorton(unsigned long* indx, unsigned long* indy, unsigned long* indz,
        unsigned long n) const;

    /** \brief Overload. */
    void deMorton(std::vector<unsigned long>* indv, unsigned long n) const;

    /** \brief Pad an integer with zeros */
    unsigned long dilate3(unsigned long x) const;

    /** \brief Remove padding */
    unsigned long contract3(unsigned long x) const;

    shared_ptr<const Bempp::Grid> m_grid;
    unsigned int m_levels;
    std::vector<std::vector<OctreeNode<CoordinateType> > > m_OctreeNodes;

    BoundingBox<CoordinateType> m_boundingBox;
};
}

#include "octree_impl.hpp"

#endif
