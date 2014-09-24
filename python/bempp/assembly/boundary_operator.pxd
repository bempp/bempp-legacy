from libcpp cimport bool as cbool


cdef extern from "bempp/assembly/boundary_operator.hpp":
    cdef cppclass c_BoundaryOperator "Bempp::BoundaryOperator" [BASIS, RESULT]:
        c_BoundaryOperator()


cdef class BoundaryOperator:
    cdef:
        void * memory
        int result_type
        int basis_type
        cbool constructed
