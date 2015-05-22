<%
from data_types import dtypes, scalar_cython_type
%>

from bempp.utils cimport Matrix
from bempp.utils.enum_types cimport TranspositionMode
from bempp.utils cimport shared_ptr
from bempp.utils cimport catch_exception
from bempp.utils cimport complex_float,complex_double
from bempp.hmat.hmatrix cimport c_HMatrix
cimport numpy as np


cdef extern from "bempp/assembly/discrete_boundary_operator.hpp" namespace "Bempp":
    cdef cppclass c_DiscreteBoundaryOperator "Bempp::DiscreteBoundaryOperator"[ValueType]:
        
        object apply(const TranspositionMode trans, object x_in) except+catch_exception

        Matrix[ValueType] asMatrix() const
        unsigned int rowCount() const
        unsigned int columnCount() const

cdef extern from "bempp/assembly/py_discrete_operator_support.hpp" namespace "Bempp":
    cdef object py_array_from_dense_operator[VALUE](const shared_ptr[const c_DiscreteBoundaryOperator[VALUE]]&)
    cdef shared_ptr[const c_HMatrix[VALUE]] py_hmatrix_from_discrete_operator[VALUE](const shared_ptr[const c_DiscreteBoundaryOperator[VALUE]]&)

cdef class DiscreteBoundaryOperatorBase:
    cdef object _dtype

cdef class DiscreteBoundaryOperator(DiscreteBoundaryOperatorBase):

% for pybasis, cybasis in dtypes.items():
    cdef shared_ptr[const c_DiscreteBoundaryOperator[${cybasis}]] _impl_${pybasis}_
% endfor

% for pyvalue, cyvalue in dtypes.items():
    cdef np.ndarray _as_matrix_${pyvalue}(self)
% endfor

cdef class SparseDiscreteBoundaryOperator(DiscreteBoundaryOperatorBase):
    cdef object _op

cdef class DenseDiscreteBoundaryOperator(DiscreteBoundaryOperator):
    cdef object _array_view
    
    cdef object _init_array_view(self)

cdef class HMatDiscreteBoundaryOperator(DiscreteBoundaryOperator):
    pass

