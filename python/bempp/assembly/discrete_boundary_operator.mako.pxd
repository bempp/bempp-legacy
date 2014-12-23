<%
from data_types import dtypes, scalar_cython_type
%>

from bempp.utils.armadillo cimport Mat
from bempp.utils.enum_types cimport TranspositionMode
from bempp.utils cimport shared_ptr
from bempp.utils cimport complex_float,complex_double
cimport numpy as np

cdef extern from "bempp/assembly/discrete_boundary_operator.hpp" namespace "Bempp":
    cdef cppclass c_DiscreteBoundaryOperator "Bempp::DiscreteBoundaryOperator"[ValueType]:
        
        void apply(const TranspositionMode trans, const Mat[ValueType]& x_in,
                Mat[ValueType]& y_inout, const ValueType alpha,
                const ValueType beta) const

        Mat[ValueType] asMatrix() const
        unsigned int rowCount() const
        unsigned int columnCount() const

cdef class DiscreteBoundaryOperatorBase:
    cdef object _value_type

% for pyvalue,cyvalue in dtypes.items():
    cdef void _apply_${pyvalue}(self,
            TranspositionMode trans, const Mat[${cyvalue}]& x_in,
            Mat[${cyvalue}]& y_inout, 
            ${scalar_cython_type(cyvalue)} alpha,
            ${scalar_cython_type(cyvalue)} beta)
% endfor

    cpdef object apply(self,np.ndarray x,np.ndarray y,object transpose,object alpha, object beta)
    

cdef class DiscreteBoundaryOperator(DiscreteBoundaryOperatorBase):

% for pybasis, cybasis in dtypes.items():
    cdef shared_ptr[const c_DiscreteBoundaryOperator[${cybasis}]] _impl_${pybasis}_
% endfor

% for pyvalue,cyvalue in dtypes.items():
    cdef void _apply_${pyvalue}(self,
            TranspositionMode trans, const Mat[${cyvalue}]& x_in,
            Mat[${cyvalue}]& y_inout, 
            ${scalar_cython_type(cyvalue)} alpha,
            ${scalar_cython_type(cyvalue)} beta)
% endfor
    
