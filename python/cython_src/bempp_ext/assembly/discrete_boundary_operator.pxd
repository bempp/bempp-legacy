from bempp_ext.utils cimport Matrix
from bempp_ext.utils.enum_types cimport TranspositionMode
from bempp_ext.utils cimport shared_ptr
from bempp_ext.utils cimport catch_exception
from bempp_ext.utils cimport complex_double

cdef extern from "bempp/assembly/discrete_boundary_operator.hpp" namespace "Bempp":
    cdef cppclass c_DiscreteBoundaryOperator "Bempp::DiscreteBoundaryOperator"[ValueType]:
        
        object apply(const TranspositionMode trans, object x_in) except+catch_exception

        Matrix[ValueType] asMatrix() const
        unsigned int rowCount() const
        unsigned int columnCount() const

cdef class DiscreteBoundaryOperatorRealExt:
    cdef shared_ptr[c_DiscreteBoundaryOperator[double]] impl_
    cdef TranspositionMode transpose_mode

cdef class DiscreteBoundaryOperatorComplexExt:
    cdef shared_ptr[c_DiscreteBoundaryOperator[complex_double]] impl_
    cdef TranspositionMode transpose_mode
