from libcpp cimport bool as cbool
from bempp.utils.armadillo cimport Mat

cdef extern from "bempp/grid/geometry.hpp" namespace "Bempp":
    cdef cppclass c_Geometry "Bempp::Geometry":
        int dim() const
        int dimWorld() const
        cbool affine() const
        int cornerCount() const
        void getCorners(Mat[double]& c) const

cdef class Geometry:
    cdef const c_Geometry * impl_
