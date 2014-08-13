from libcpp cimport bool as cbool

cdef extern from "bempp/fiber/quadrature_options.hpp" namespace "Fiber":
    cdef cppclass c_QuadratureOptions "Fiber::QuadratureOptions":
        c_QuadratureOptions()
        c_QuadratureOptions(int order, cbool relative)
        void setAbsoluteQuadratureOrder(int order)
        void setRelativeQuadratureOrder(int offset)
        int quadratureOrder(int defaultOrder) const

cdef class QuadratureOptions:
    cdef:
        c_QuadratureOptions impl
        cbool __is_relative
        cbool __is_frozen
