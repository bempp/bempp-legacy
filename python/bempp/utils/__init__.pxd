from bempp.utils.shared_ptr cimport shared_ptr
from bempp.utils.unique_ptr cimport unique_ptr

cdef extern from "bempp/utils/utils.h" namespace "Bempp":
    void catch_exception()
