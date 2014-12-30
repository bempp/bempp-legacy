<%
from data_types import dtypes,scalar_cython_type
%>
cimport cython
cimport numpy as np
import numpy as np
from bempp.utils cimport complex_float,complex_double
from cython.operator cimport dereference as deref

% for pyvalue,cyvalue in dtypes.items():
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray armadillo_to_np_${pyvalue}(const Mat[${cyvalue}]& x):
    
    cdef int rows = x.n_rows
    cdef int cols = x.n_cols
    cdef int j
    cdef int i

    cdef np.ndarray[${scalar_cython_type(cyvalue)},ndim=2,mode='fortran'] res = np.empty((rows,cols),dtype="${pyvalue}",order='F')

    for j in range(cols):
        for i in range(rows):
            res[i,j] = deref(<${scalar_cython_type(cyvalue)}*>&x.at(i,j))
    return res
% endfor


% for pyvalue,cyvalue in dtypes.items():
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray armadillo_col_to_np_${pyvalue}(const Col[${cyvalue}]& x):
    
    cdef int rows = x.n_rows
    cdef int i

    cdef np.ndarray[${scalar_cython_type(cyvalue)},ndim=1,mode='fortran'] res = np.empty(rows,dtype="${pyvalue}",order='F')

    for i in range(rows):
        res[i] = deref(<${scalar_cython_type(cyvalue)}*>&x.at(i))
    return res
% endfor
