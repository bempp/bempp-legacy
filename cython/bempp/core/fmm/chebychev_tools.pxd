from bempp.core.utils.shared_ptr cimport shared_ptr
from bempp.core.utils.eigen cimport Vector, Matrix

from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "bempp/fmm/chebychev_tools.hpp":
    cdef cppclass c_ChebychevTools "Fmm::ChebychevTools":
        c_ChebychevTools(int order)
        const Vector[double]& chebychevNodes() const
        void evaluateInterpolationPolynomial(const Vector[double]& weights, const Vector[double]& evaluationPoints,
                                             Vector[double] result) const
        Matrix[double] childInterpolationMatrix(double ratio) const
        Vector[double] derivativeWeights(const Vector[double]& weights) const
        Vector[double] derivativeWeights3d(const Vector[double]& weights, int direction) const


cdef class ChebychevTools:
    cdef shared_ptr[c_ChebychevTools] impl_

