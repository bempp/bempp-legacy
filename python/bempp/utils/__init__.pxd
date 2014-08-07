from bempp.utils.shared_ptr cimport shared_ptr

cdef extern from "bempp/utils/utils.h" namespace "Bempp":
    void catch_exception()
