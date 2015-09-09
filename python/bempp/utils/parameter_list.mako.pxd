from bempp.utils cimport shared_ptr
from libcpp.string cimport string 
from libcpp cimport bool as cbool
from libcpp.vector cimport vector


cdef extern from "bempp/common/types.hpp" namespace "Bempp":

    cdef cppclass c_ParameterList "Bempp::ParameterList":
        c_ParameterList()
        c_ParameterList& assign "operator="(const c_ParameterList&)
        void put_bool "put<bool>" (char*, cbool)
        void put_string "put<std::string>" (char*, string)
        void put_int "put<int>" (char*, int)
        void put_double "put<double>" (char*, double)
        string get_string "get<std::string>" (char*)
        int get_int "get<int>" (char*)
        double get_double "get<double>" (char*)
        cbool get_bool "get<bool>" (char*)

cdef extern from "<sstream>" namespace "std":

    cdef cppclass ostringstream:
        ostringstream() except +
        string str()
        void str(string)
        void clear()

cdef extern from "bempp/common/common.hpp" namespace "boost::property_tree::json_parser":

    void write_file "write_json" (string, c_ParameterList)
    void write_stream "write_json"(ostringstream, c_ParameterList, cbool)



cdef class _NearField:
    cdef c_ParameterList* impl_
    cdef _QuadratureParameterList base

cdef class _MediumField:
    cdef c_ParameterList* impl_
    cdef _QuadratureParameterList base

cdef class _FarField:
    cdef c_ParameterList* impl_
    cdef _QuadratureParameterList base

cdef class _AssemblyParameterList:
    cdef c_ParameterList* impl_
    cdef ParameterList base

cdef class _QuadratureParameterList:
    cdef c_ParameterList* impl_
    cdef _NearField  _near
    cdef _MediumField _medium
    cdef _FarField _far
    cdef ParameterList base

cdef class _HMatParameterList:
    cdef c_ParameterList* impl_
    cdef ParameterList base


cdef class ParameterList:
    cdef c_ParameterList* impl_
    cdef ostringstream* outputter_
    cdef _AssemblyParameterList _assembly
    cdef _QuadratureParameterList _quadrature
    cdef _HMatParameterList _hmat
