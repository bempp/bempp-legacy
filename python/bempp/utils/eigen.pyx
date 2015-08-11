cimport cython
cimport numpy as np
import numpy as np
from bempp.utils cimport complex_float,complex_double
from cython.operator cimport dereference as deref

cdef Matrix[double] np_to_eigen_matrix_float64(np.ndarray x):
    cdef double[::1,:] buf = np.require(x,dtype='float64',requirements=['A','F'])
    return copy_buf_to_mat[double](<double*>&buf[0,0],buf.shape[0],buf.shape[1])

cdef Matrix[complex_double] np_to_eigen_matrix_complex128(np.ndarray x):
    cdef double complex[::1,:] buf = np.require(x,dtype='complex128',requirements=['A','F'])
    return copy_buf_to_mat[complex_double](<complex_double*>&buf[0,0],buf.shape[0],buf.shape[1])

cdef Vector[double] np_to_eigen_vector_float64(np.ndarray x):
    cdef double[::1] buf = np.require(x,dtype='float64',requirements=['A','F'])
    return copy_buf_to_vec[double](<double*>&buf[0],buf.shape[0])

cdef Vector[complex_double] np_to_eigen_vector_complex128(np.ndarray x):
    cdef double complex[::1] buf = np.require(x,dtype='complex128',requirements=['A','F'])
    return copy_buf_to_vec[complex_double](<complex_double*>&buf[0],buf.shape[0])

@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray eigen_matrix_to_np_float64(const Matrix[double]& x):
    
    cdef int rows = x.rows()
    cdef int cols = x.cols()
    cdef int j
    cdef int i

    cdef np.ndarray[double,ndim=2,mode='fortran'] res = np.empty((rows,cols),dtype="float64",order='F')

    for j in range(cols):
        for i in range(rows):
            res[i,j] = deref(<double*>&x.value(i,j))
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray eigen_matrix_to_np_complex128(const Matrix[complex_double]& x):
    
    cdef int rows = x.rows()
    cdef int cols = x.cols()
    cdef int j
    cdef int i

    cdef np.ndarray[double complex,ndim=2,mode='fortran'] res = np.empty((rows,cols),dtype="complex128",order='F')

    for j in range(cols):
        for i in range(rows):
            res[i,j] = deref(<double complex*>&x.value(i,j))
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray eigen_vector_to_np_float64(const Vector[double]& x):
    
    cdef int rows = x.rows()
    cdef int i

    cdef np.ndarray[double,ndim=1,mode='fortran'] res = np.empty(rows,dtype="float64",order='F')

    for i in range(rows):
        res[i] = deref(<double*>&x.value(i))
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray eigen_vector_to_np_complex128(const Vector[complex_double]& x):
    
    cdef int rows = x.rows()
    cdef int i

    cdef np.ndarray[double complex,ndim=1,mode='fortran'] res = np.empty(rows,dtype="complex128",order='F')

    for i in range(rows):
        res[i] = deref(<double complex*>&x.value(i))
    return res

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

def _test_eigen_matrix_conversion_float64(x):
    return eigen_matrix_to_np_float64(np_to_eigen_matrix_float64(x))

def _test_eigen_vector_conversion_float64(x):
    return eigen_vector_to_np_float64(np_to_eigen_vector_float64(x))

def _test_eigen_matrix_conversion_complex128(x):
    return eigen_matrix_to_np_complex128(np_to_eigen_matrix_complex128(x))

def _test_eigen_vector_conversion_complex128(x):
    return eigen_vector_to_np_complex128(np_to_eigen_vector_complex128(x))
