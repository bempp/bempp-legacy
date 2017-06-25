from bempp.core.utils.shared_ptr cimport shared_ptr
from bempp.core.grid.grid cimport c_Grid
from bempp.core.common.bounding_box cimport BoundingBox


cdef extern from "bempp/fmm/octree_node.hpp":
    cdef cppclass c_OctreeNode "Fmm::OctreeNode<double>":
        pass

cdef extern from "bempp/fmm/octree.hpp":
    cdef cppclass c_Octree "Fmm::Octree<double>":
        c_Octree(shared_ptr[c_Grid], unsigned int levels)
        BoundingBox getBoundingBox()
        unsigned long getParent(unsigned long n)
        unsigned long getFirstChild(unsigned long n)
        unsigned long getLastChild(unsigned long n)
        unsigned long getNodesPerSide(unsigned long level)
        unsigned long getNodesPerLevel(unsigned long level)

cdef class OctreeNode:
    cdef shared_ptr[c_OctreeNode] impl_

cdef class Octree:
    cdef shared_ptr[c_Octree] impl_
