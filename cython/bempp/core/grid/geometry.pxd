from libcpp cimport bool as cbool
from bempp.core.utils cimport Matrix, RowVector

cdef extern from "bempp/grid/geometry.hpp" namespace "Bempp":
    cdef cppclass c_Geometry "Bempp::Geometry":
        int dim() const
        int dimWorld() const
        cbool affine() const
        int cornerCount() const
        void getCorners(Matrix[double]& c) const
        double volume() const
        void getIntegrationElements(const Matrix[double]& local,
                RowVector[double]& int_element)

from bempp.core.grid.entity cimport Entity0
from bempp.core.grid.entity cimport Entity1
from bempp.core.grid.entity cimport Entity2

cdef class Geometry0:
    cdef const c_Geometry * impl_
    cdef Entity0 _entity

cdef class Geometry1:
    cdef const c_Geometry * impl_
    cdef Entity1 _entity

cdef class Geometry2:
    cdef const c_Geometry * impl_
    cdef Entity2 _entity
