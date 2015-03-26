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
cdef np.ndarray eigen_matrix_to_np_${pyvalue}(const Matrix[${cyvalue}]& x):
    
    cdef int rows = x.rows()
    cdef int cols = x.cols()
    cdef int j
    cdef int i

    cdef np.ndarray[${scalar_cython_type(cyvalue)},ndim=2,mode='fortran'] res = np.empty((rows,cols),dtype="${pyvalue}",order='F')

    for j in range(cols):
        for i in range(rows):
            res[i,j] = deref(<${scalar_cython_type(cyvalue)}*>&x(i,j))
    return res
% endfor


% for pyvalue,cyvalue in dtypes.items():
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray eigen_vector_to_np_${pyvalue}(const Vector[${cyvalue}]& x):
    
    cdef int rows = x.rows()
    cdef int i

    cdef np.ndarray[${scalar_cython_type(cyvalue)},ndim=1,mode='fortran'] res = np.empty(rows,dtype="${pyvalue}",order='F')

    for i in range(rows):
        res[i] = deref(<${scalar_cython_type(cyvalue)}*>&x(i))
    return res
% endfor


@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray eigen_to_np_int(const Matrix[int]& x):
    
    cdef int rows = x.rows()
    cdef int cols = x.cols()
    cdef int j
    cdef int i

    cdef np.ndarray[int,ndim=2,mode='fortran'] res = np.empty((rows,cols),dtype="intc",order='F')

    for j in range(cols):
        for i in range(rows):
            res[i,j] = x(i,j)
    return res
