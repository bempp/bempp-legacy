from bempp.utils.shared_ptr cimport shared_ptr
from bempp.utils.unique_ptr cimport unique_ptr
from bempp.utils.parameter_list cimport ParameterList

cdef extern from "bempp/utils/utils.hpp" namespace "Bempp":
    void catch_exception()
