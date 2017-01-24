cdef extern from "structure.h":
    ctypedef struct Structure:
        pass

    void init_structure(Structure *);
    void dealloc_structure(Structure *);

cdef class Structure:
    cdef Structure* __internal_structure
