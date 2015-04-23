<%
from data_types import dtypes,scalar_cython_type
%>

from bempp.utils.shared_ptr cimport shared_ptr, static_pointer_cast
from bempp.utils.unique_ptr cimport unique_ptr

from bempp.utils.parameter_list cimport ParameterList

cdef extern from "bempp/utils/py_utils.hpp" namespace "Bempp":
    void catch_exception()

cdef extern from "<complex>" namespace "std":
    cdef cppclass cpp_complex "std::complex"[T]:
        cpp_complex()
        cpp_complex(T alpha,T beta)

ctypedef cpp_complex[float] complex_float
ctypedef cpp_complex[double] complex_double

from bempp.utils.eigen cimport Vector, Matrix

