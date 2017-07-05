from cython.operator cimport dereference as deref
from cython.operator cimport address
from bempp.core.utils.shared_ptr cimport shared_ptr
from bempp.core.grid.grid cimport Grid, c_Grid
from bempp.core.utils cimport eigen_vector_to_np_float64

cdef class Octree:
    """Define an Octree over a given grid."""

    def __cinit__(self, Grid grid, int levels):
        pass

    def __init__(self, Grid grid, int levels):
        self.impl_.assign(shared_ptr[c_Octree](new c_Octree(grid.impl_, levels)))

    def __dealloc__(self):
        self.impl_.reset()

    def parent(self, n):
        """Return parent index."""
        return deref(self.impl_).getParent(n)

    def children(self, n):
        """Return [first_child, last_child+1]."""
        return (deref(self.impl_).getFirstChild(n), deref(self.impl_).getLastChild(n) + 1)

    def nodes_per_side(self, level):
        """Get the number of nodes on each side in this level."""
        return deref(self.impl_).getNodesPerSide(level)

    def nodes_per_level(self, level):
        """ Get the total number of nodes in the level."""
        return deref(self.impl_).getNodesPerLevel(level)

    def cube_width(self, level):
        """ Get the width of a cube on a given level."""
        return deref(self.impl_).cubeWidth(level)

    def extended_cube_width(self, level):
        """ Get the extended width of a cube on a given level."""
        return deref(self.impl_).extendedCubeWidth(level)

    def cube_bounds(self, unsigned long node_index, unsigned int level):
        """Return lower and upper bounds for a given cube."""

        cdef Vector[double] lbound
        cdef Vector[double] ubound

        deref(self.impl_).cubeBounds(node_index, level, lbound, ubound)
        return (eigen_vector_to_np_float64(lbound),
                eigen_vector_to_np_float64(ubound))

    def extended_cube_bounds(self, unsigned long node_index, unsigned int level):
        """Return lower and upper bounds for a given cube."""

        cdef Vector[double] lbound
        cdef Vector[double] ubound

        deref(self.impl_).extendedCubeBounds(node_index, level, lbound, ubound)
        return (eigen_vector_to_np_float64(lbound),
                eigen_vector_to_np_float64(ubound))
        
    def is_empty(self, unsigned long node_index, unsigned int level):
        """Return if a given node is empty."""
        return deref(self.impl_).isEmpty(node_index, level)

    def leaf_cube_entities(self, unsigned long node_index):
        """Return all element indices associated with a leaf cube."""

        return deref(self.impl_).getLeafCubeEntities(node_index)


    property bounding_box:

        def __get__(self):
            cdef BoundingBox bbox = deref(self.impl_).getBoundingBox()
            res = {'lbound':[bbox.lbound.x, bbox.lbound.y, bbox.lbound.z],
                    'ubound':[bbox.ubound.x, bbox.ubound.y, bbox.ubound.z]}
            return res

    property levels:
        """Return the number of levels."""

        def __get__(self):
            return deref(self.impl_).levels()
