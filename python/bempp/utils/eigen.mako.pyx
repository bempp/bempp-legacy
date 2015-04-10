<%
from data_types import dtypes,scalar_cython_type
%>
cimport cython
cimport numpy as np
import numpy as np
from bempp.utils cimport complex_float,complex_double
from cython.operator cimport dereference as deref

% for pyvalue,cyvalue in dtypes.items():

cdef Matrix[${cyvalue}] np_to_eigen_matrix_${pyvalue}(np.ndarray x):

    cdef ${scalar_cython_type(cyvalue)}[::1,:] buf = np.require(x,dtype='${pyvalue}',requirements=['A','F'])

    return copy_buf_to_mat[${cyvalue}](<${cyvalue}*>&buf[0,0],buf.shape[0],buf.shape[1])

cdef Vector[${cyvalue}] np_to_eigen_vector_${pyvalue}(np.ndarray x):
    cdef ${scalar_cython_type(cyvalue)}[::1] buf = np.require(x,dtype='${pyvalue}',requirements=['A','F'])

    return copy_buf_to_vec[${cyvalue}](<${cyvalue}*>&buf[0],buf.shape[0])


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
            res[i,j] = deref(<${scalar_cython_type(cyvalue)}*>&x.value(i,j))
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
        res[i] = deref(<${scalar_cython_type(cyvalue)}*>&x.value(i))
    return res
% endfor


@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray eigen_vector_to_np_int(const Vector[int]& x):
    
    cdef int rows = x.rows()
    cdef int i

    cdef np.ndarray[int,ndim=1,mode='fortran'] res = np.empty(rows,dtype="intc",order='F')

    for i in range(rows):
        res[i] = x.value(i)
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray eigen_matrix_to_np_int(const Matrix[int]& x):
    
    cdef int rows = x.rows()
    cdef int cols = x.cols()
    cdef int j
    cdef int i

    cdef np.ndarray[int,ndim=2,mode='fortran'] res = np.empty((rows,cols),dtype="intc",order='F')

    for j in range(cols):
        for i in range(rows):
            res[i,j] = x.value(i,j)
    return res

% for pyvalue,cyvalue in dtypes.items():

def _test_eigen_matrix_conversion_${pyvalue}(x):

    return eigen_matrix_to_np_${pyvalue}(np_to_eigen_matrix_${pyvalue}(x))

def _test_eigen_vector_conversion_${pyvalue}(x):

    return eigen_vector_to_np_${pyvalue}(np_to_eigen_vector_${pyvalue}(x))

% endfor
