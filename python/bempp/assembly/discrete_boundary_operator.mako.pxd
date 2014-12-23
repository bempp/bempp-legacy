<%
from data_types import dtypes
%>

from bempp.utils.armadillo cimport Mat
from bempp.utils.enum_types cimport TranspositionMode
from bempp.utils cimport shared_ptr
from bempp.utils cimport complex_float,complex_double


cdef extern from "bempp/assembly/discrete_boundary_operator.hpp" namespace "Bempp":
    cdef cppclass c_DiscreteBoundaryOperator "Bempp::DiscreteBoundaryOperator"[ValueType]:
        
        void apply(const TranspositionMode trans, const Mat[ValueType]& x_in,
                Mat[ValueType]& y_inout, const ValueType alpha,
                const ValueType beta) const

        Mat[ValueType] asMatrix() const
        unsigned int rowCount() const
        unsigned int columnCount() const

cdef class DiscreteBoundaryOperator:
    cdef object _value_type

% for pybasis, cybasis in dtypes.items():
    cdef shared_ptr[const c_DiscreteBoundaryOperator[${cybasis}]] _impl_${pybasis}_
% endfor
    
