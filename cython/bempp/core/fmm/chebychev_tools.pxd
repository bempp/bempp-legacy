from bempp.core.utils.shared_ptr cimport shared_ptr
from bempp.core.utils.eigen cimport Vector

from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "bempp/fmm/chebychev_tools.hpp":
    cdef cppclass c_ChebychevTools "Fmm::ChebychevTools":
        c_ChebychevTools(int order)
        void chebychevNodes(const Vector[double]& nodes)


cdef class ChebychevTools:
    cdef shared_ptr[c_ChebychevTools] impl_

