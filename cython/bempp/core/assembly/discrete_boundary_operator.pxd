from bempp.core.utils cimport Matrix
from bempp.core.utils.enum_types cimport TranspositionMode
from bempp.core.utils cimport shared_ptr
from bempp.core.utils cimport catch_exception
from bempp.core.utils cimport complex_double

cdef extern from "bempp/assembly/discrete_boundary_operator.hpp" namespace "Bempp":
    cdef cppclass c_DiscreteBoundaryOperator "Bempp::DiscreteBoundaryOperator"[ValueType]:
        
        Matrix[ValueType] apply(const TranspositionMode trans, Matrix[ValueType] x_in) except+catch_exception

        Matrix[ValueType] asMatrix() const
        unsigned int rowCount() const
        unsigned int columnCount() const

cdef class RealDiscreteBoundaryOperator:
    cdef shared_ptr[const c_DiscreteBoundaryOperator[double]] impl_
    cdef TranspositionMode transpose_mode

cdef class ComplexDiscreteBoundaryOperator:
    cdef shared_ptr[const c_DiscreteBoundaryOperator[complex_double]] impl_
    cdef TranspositionMode transpose_mode

cdef extern from "bempp/core/assembly/discrete_operator_conversion.hpp" namespace "Bempp":
    cdef object c_convert_to_sparse "Bempp::py_get_sparse_from_discrete_operator<double>"(const shared_ptr[c_DiscreteBoundaryOperator[double]]&)
    cdef Matrix[double] c_convert_to_dense_double "Bempp::eigen_matrix_from_dense_operator<double>"(const shared_ptr[c_DiscreteBoundaryOperator[double]]&)
    cdef Matrix[complex_double] c_convert_to_dense_complex "Bempp::eigen_matrix_from_dense_operator<std::complex<double> >"(const shared_ptr[c_DiscreteBoundaryOperator[complex_double]]&)




