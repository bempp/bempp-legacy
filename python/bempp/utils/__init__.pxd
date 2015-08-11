from bempp.utils.shared_ptr cimport shared_ptr, static_pointer_cast, const_pointer_cast, reverse_const_pointer_cast
from bempp.utils.unique_ptr cimport unique_ptr

from bempp.utils.parameter_list cimport ParameterList

cdef extern from "bempp/utils/py_utils.hpp" namespace "Bempp":
    void catch_exception()

cdef extern from "<complex>" namespace "std":
    cdef cppclass cpp_complex "std::complex"[T]:
        cpp_complex()
        cpp_complex(T alpha,T beta)

ctypedef cpp_complex[double] complex_double

from bempp.utils.eigen cimport Vector, Matrix
from bempp.utils.eigen cimport eigen_matrix_to_np_float64
from bempp.utils.eigen cimport eigen_matrix_to_np_complex128
from bempp.utils.eigen cimport eigen_vector_to_np_float64
from bempp.utils.eigen cimport eigen_vector_to_np_complex128
from bempp.utils.eigen cimport np_to_eigen_matrix_float64
from bempp.utils.eigen cimport np_to_eigen_matrix_complex128
from bempp.utils.eigen cimport np_to_eigen_vector_float64
from bempp.utils.eigen cimport np_to_eigen_vector_complex128
from bempp.utils.eigen cimport eigen_matrix_to_np_int
from bempp.utils.eigen cimport eigen_vector_to_np_int

