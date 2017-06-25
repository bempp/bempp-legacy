cdef extern from "bempp/common/types.hpp":
    cdef cppclass Point "Bempp::Point3D<double>":
        double x
        double y
        double z
