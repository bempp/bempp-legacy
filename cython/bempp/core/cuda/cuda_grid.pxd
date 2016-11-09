from bempp.core.utils cimport shared_ptr, unique_ptr

cdef extern from "bempp/cuda/cuda_grid.hpp" namespace "Bempp":
    cdef cppclass c_CudaGrid "Bempp::CudaGrid":
        pass


cdef class CudaGrid:
    ## Holds pointer to C++ implementation
    cdef shared_ptr[c_CudaGrid] impl_
