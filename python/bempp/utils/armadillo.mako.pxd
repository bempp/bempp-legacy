<%
from data_types import dtypes,scalar_cython_type
%>


from libcpp cimport bool as cbool
cimport numpy as np
from bempp.utils cimport complex_float,complex_double

cdef extern from "<armadillo>" namespace "arma":
    ctypedef size_t uword
    cdef cppclass Col[T]:
        Col(T*, uword, cbool copy_aux_mem, cbool strict) nogil
        Col() nogil
        T* memptr()
        T& value "operator()"(int i) # bounds checking
        T& at(int i) # No bounds checking
        int n_rows
        int n_cols

    cdef cppclass Mat[T]:
        Mat(T*, uword, uword, cbool copy_aux_mem, cbool strict) nogil
        T& at(int i, int j) # no bounds checking
        T& value "operator()"(int i, int j) # bounds checking
        Mat() nogil
        T* memptr()
        int n_rows
        int n_cols

% for pyvalue,cyvalue in dtypes.items():
cdef np.ndarray armadillo_to_np_${pyvalue}(const Mat[${cyvalue}]& x)
% endfor
