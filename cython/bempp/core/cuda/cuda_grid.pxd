from bempp.core.utils cimport shared_ptr, unique_ptr

cdef extern from "bempp/cuda/cuda_grid.hpp" namespace "Bempp":
    cdef cppclass c_CudaGridFloat "Bempp::CudaGrid<float>":
        pass

cdef extern from "bempp/cuda/cuda_grid.hpp" namespace "Bempp":
    cdef cppclass c_CudaGridDouble "Bempp::CudaGrid<double>":
        pass

cdef class CudaGridFloat:
    ## Holds pointer to C++ implementation
    cdef shared_ptr[c_CudaGridFloat] impl_

cdef class CudaGridDouble:
    ## Holds pointer to C++ implementation
    cdef shared_ptr[c_CudaGridDouble] impl_
