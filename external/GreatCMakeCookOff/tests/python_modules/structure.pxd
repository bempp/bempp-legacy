cdef extern from "structure.h":
    ctypedef struct CStructure "Structure":
        int meaning_of_life
        const char * message

    void init_structure(CStructure *_self)
    void dealloc_structure(CStructure *_self)

cdef class Structure:
    cdef CStructure* cdata
