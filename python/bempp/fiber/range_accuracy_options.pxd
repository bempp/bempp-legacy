from libcpp cimport bool as cbool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from bempp.fiber.quadrature_options cimport c_QuadratureOptions


cdef class RangeAccuracyOptions(dict):
    ## Accuracy of floating points indices
    cdef public double __tolerance__
    ## Wether the object can be modified or not
    cdef cbool __is_frozen
    cdef double __index(self, key)
    cdef toggle_freeze(self, value=?)
    cdef void to_cpp(self,
            vector[pair[double, c_QuadratureOptions]] &quadops) except *
