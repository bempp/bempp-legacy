from bempp.core.utils.shared_ptr cimport shared_ptr
from bempp.core.grid.grid cimport c_Grid
from bempp.core.common.bounding_box cimport BoundingBox
from bempp.core.utils.eigen cimport Vector


cdef extern from "bempp/fmm/octree.hpp":
    cdef cppclass c_Octree "Fmm::Octree":
        c_Octree(shared_ptr[c_Grid], unsigned int levels)
        BoundingBox getBoundingBox()
        unsigned int levels()
        unsigned long getParent(unsigned long nodeIndex)
        unsigned long getFirstChild(unsigned long nodeIndex)
        unsigned long getLastChild(unsigned long nodeIndex)
        unsigned long getNodesPerSide(unsigned int level)
        unsigned long getNodesPerLevel(unsigned int level)
        unsigned long isEmpty(unsigned long nodeIndex, unsigned int level)
        double cubeWidth(unsigned int level);
        double extendedCubeWidth(unsigned int level);
        void cubeBounds(unsigned long nodeIndex, unsigned int level, 
                Vector[double]& lbound, Vector[double]& ubound)
        void extendedCubeBounds(unsigned long nodeIndex, unsigned int level, 
                Vector[double]& lbound, Vector[double]& ubound)
        

cdef class Octree:
    cdef shared_ptr[c_Octree] impl_
