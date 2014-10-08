from libcpp cimport bool as cbool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from bempp.fiber.range_accuracy_options cimport RangeAccuracyOptions
from bempp.fiber.quadrature_options cimport QuadratureOptions, \
        c_QuadratureOptions

cdef extern from "bempp/fiber/accuracy_options.hpp" namespace "Fiber":
    cdef cppclass c_AccuracyOptions "Fiber::AccuracyOptionsEx":
        c_AccuracyOptions()
        void setDoubleRegular(const vector[pair[double, c_QuadratureOptions]]&)
        void setSingleRegular(const vector[pair[double, c_QuadratureOptions]]&)
        void setDoubleSingular(int, cbool)
        const c_QuadratureOptions& singleRegular(double) const
        const c_QuadratureOptions& doubleRegular(double) const
        const c_QuadratureOptions& doubleSingular() const

cdef class AccuracyOptions:
    cdef:
       RangeAccuracyOptions __single_regular
       RangeAccuracyOptions __double_regular
       QuadratureOptions __double_singular
    cdef toggle_freeze(self, value=?)
    cdef void to_cpp(self, c_AccuracyOptions &) except*
