from bempp.core.utils cimport shared_ptr, unique_ptr

cdef extern from "bempp/cuda/cuda_grid.hpp" namespace "Bempp":
    cdef cppclass c_CudaGridSingle "Bempp::CudaGrid<float>":
        pass

cdef extern from "bempp/cuda/cuda_grid.hpp" namespace "Bempp":
    cdef cppclass c_CudaGridDouble "Bempp::CudaGrid<double>":
        pass

cdef class CudaGridSingle:
    ## Holds pointer to C++ implementation
    cdef shared_ptr[c_CudaGridSingle] impl_

cdef class CudaGridDouble:
    ## Holds pointer to C++ implementation
    cdef shared_ptr[c_CudaGridDouble] impl_
