from bempp.core.utils.shared_ptr cimport shared_ptr
from bempp.core.fmm.octree cimport c_Octree

from libcpp cimport bool
from libcpp.vector cimport vector


cdef extern from "bempp/fmm/interaction_list.hpp":
    cdef cppclass c_InteractionList "Fmm::InteractionList":
        c_InteractionList(const c_Octree&, unsigned long, unsigned int)
        void next()
        bool finished() 
        unsigned long& operator*()

        
cdef class InteractionList:
    cdef shared_ptr[c_InteractionList] impl_
