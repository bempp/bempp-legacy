cdef class BoundaryOperator:
    cdef:
        void * memory
        int result_type
        int basis_type
