from bempp.utils cimport shared_ptr
from libcpp.string cimport string 
from libcpp cimport bool as cbool
from libcpp.vector cimport vector


cdef extern from "bempp/common/types.hpp" namespace "Bempp":

    cdef cppclass c_ParameterList "Bempp::ParameterList":
        c_ParameterList()
        c_ParameterList& assign "operator="(const c_ParameterList&)


cdef class ParameterList:
    cdef c_ParameterList* impl_
