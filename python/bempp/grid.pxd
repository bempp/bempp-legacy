from bempp.utils cimport shared_ptr

cdef extern from "bempp/grid/grid.hpp" namespace "Bempp":
    cdef cppclass c_Grid "Bempp::Grid":
        int dim() const
        int dimWorld() const
        int maxLevel() const
        int topology() const

cdef class Grid:
    ## Holds pointer to C++ implementation
    cdef shared_ptr[c_Grid] impl_
