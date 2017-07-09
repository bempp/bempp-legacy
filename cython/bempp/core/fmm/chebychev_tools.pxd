from bempp.core.utils.shared_ptr cimport shared_ptr
from bempp.core.utils.eigen cimport Vector, Matrix

from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "bempp/fmm/chebychev_tools.hpp":
    cdef cppclass c_ChebychevTools "Fmm::ChebychevTools":
        c_ChebychevTools(int order)
        const Vector[double]& chebychevNodes() const
        const Matrix[double]& chebychevPolValuesAtNodes() const

cdef class ChebychevTools:
    cdef shared_ptr[c_ChebychevTools] impl_

