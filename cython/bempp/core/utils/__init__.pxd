from bempp.core.utils.shared_ptr cimport shared_ptr, static_pointer_cast, const_pointer_cast, reverse_const_pointer_cast
from bempp.core.utils.unique_ptr cimport unique_ptr
from bempp.core.utils.complex cimport complex_double, cpp_complex
from bempp.core.utils.signal_slot_interface cimport Connection, SlotInterface

from bempp.core.utils.parameter_list cimport ParameterList, c_ParameterList

cdef extern from "bempp/core/utils/py_utils.hpp" namespace "Bempp":
    void catch_exception()


from bempp.core.utils.eigen cimport Vector, Matrix
from bempp.core.utils.eigen cimport eigen_matrix_to_np_float64
from bempp.core.utils.eigen cimport eigen_matrix_to_np_complex128
from bempp.core.utils.eigen cimport eigen_vector_to_np_float64
from bempp.core.utils.eigen cimport eigen_vector_to_np_complex128
from bempp.core.utils.eigen cimport np_to_eigen_matrix_float64
from bempp.core.utils.eigen cimport np_to_eigen_matrix_complex128
from bempp.core.utils.eigen cimport np_to_eigen_vector_float64
from bempp.core.utils.eigen cimport np_to_eigen_vector_complex128
from bempp.core.utils.eigen cimport eigen_matrix_to_np_int
from bempp.core.utils.eigen cimport eigen_vector_to_np_int

