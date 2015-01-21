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

cdef extern from "bempp/assembly/py_discrete_operator_support.hpp" namespace "Bempp":
    cdef object py_array_from_dense_operator[VALUE](const shared_ptr[const c_DiscreteBoundaryOperator[VALUE]]&)

cdef class DiscreteBoundaryOperatorBase:
    cdef object _value_type

% for pyvalue,cyvalue in dtypes.items():
    cdef void _apply_${pyvalue}(self,
            TranspositionMode trans, 
            np.ndarray[${scalar_cython_type(cyvalue)},ndim=2,mode='fortran'] x_in, 
            np.ndarray[${scalar_cython_type(cyvalue)},ndim=2,mode='fortran'] y_inout, 
            ${scalar_cython_type(cyvalue)} alpha,
            ${scalar_cython_type(cyvalue)} beta)
% endfor


cdef class DiscreteBoundaryOperator(DiscreteBoundaryOperatorBase):

% for pybasis, cybasis in dtypes.items():
    cdef shared_ptr[const c_DiscreteBoundaryOperator[${cybasis}]] _impl_${pybasis}_
% endfor

% for pyvalue,cyvalue in dtypes.items():
    cdef void _apply_${pyvalue}(self,
            TranspositionMode trans, 
            np.ndarray[${scalar_cython_type(cyvalue)},ndim=2,mode='fortran'] x_in, 
            np.ndarray[${scalar_cython_type(cyvalue)},ndim=2,mode='fortran'] y_inout, 
            ${scalar_cython_type(cyvalue)} alpha,
            ${scalar_cython_type(cyvalue)} beta)
    cdef np.ndarray _as_matrix_${pyvalue}(self)
% endfor

cdef class SparseDiscreteBoundaryOperator(DiscreteBoundaryOperatorBase):
    cdef object _op
    cdef object _op_transpose

cdef class DenseDiscreteBoundaryOperator(DiscreteBoundaryOperator):
    cdef object _array_view
    
    cdef object _init_array_view(self)
