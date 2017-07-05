from bempp.core.common.point cimport Point


cdef extern from "bempp/common/bounding_box.hpp":
    cdef cppclass BoundingBox "Bempp::BoundingBox<double>":
        Point lbound
        Point ubound
        Point reference


