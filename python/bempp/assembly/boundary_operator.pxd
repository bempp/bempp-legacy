from libcpp cimport bool as cbool

cdef class BoundaryOperator:
    cdef:
        void * memory
        int result_type
        int basis_type
        cbool constructed
