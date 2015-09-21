from bempp_ext.utils.shared_ptr cimport shared_ptr, static_pointer_cast, const_pointer_cast, reverse_const_pointer_cast
from bempp_ext.utils.unique_ptr cimport unique_ptr
from bempp_ext.utils.complex cimport complex_double, cpp_complex
from bempp_ext.utils.signal_slot_interface cimport Connection, SlotInterface

from bempp_ext.utils.parameter_list cimport ParameterList, c_ParameterList

cdef extern from "bempp_ext/utils/py_utils.hpp" namespace "Bempp":
    void catch_exception()


from bempp_ext.utils.eigen cimport Vector, Matrix
from bempp_ext.utils.eigen cimport eigen_matrix_to_np_float64
from bempp_ext.utils.eigen cimport eigen_matrix_to_np_complex128
from bempp_ext.utils.eigen cimport eigen_vector_to_np_float64
from bempp_ext.utils.eigen cimport eigen_vector_to_np_complex128
from bempp_ext.utils.eigen cimport np_to_eigen_matrix_float64
from bempp_ext.utils.eigen cimport np_to_eigen_matrix_complex128
from bempp_ext.utils.eigen cimport np_to_eigen_vector_float64
from bempp_ext.utils.eigen cimport np_to_eigen_vector_complex128
from bempp_ext.utils.eigen cimport eigen_matrix_to_np_int
from bempp_ext.utils.eigen cimport eigen_vector_to_np_int

