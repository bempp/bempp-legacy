cimport cython
cimport numpy as np
import numpy as np


@cython.boundscheck(False)
@cython.wraparound(False)
cdef object _3d_array_to_numpy(_3dArray array):

    cdef size_t n0 = array.extent(0)
    cdef size_t n1 = array.extent(1)
    cdef size_t n2 = array.extent(2)

    cdef int n = n0 * n1 * n2 
    cdef int i

    cdef double* data = array.begin()

    cdef np.ndarray[double,ndim=1,mode='fortran'] res = np.empty(n,dtype="float64",order='F')

    for i in range(n):
        res[i] = data[i]

    return res.reshape((n0, n1, n2), order='F')

@cython.boundscheck(False)
@cython.wraparound(False)
cdef object _4d_array_to_numpy(_4dArray array):

    cdef size_t n0 = array.extent(0)
    cdef size_t n1 = array.extent(1)
    cdef size_t n2 = array.extent(2)
    cdef size_t n3 = array.extent(3)

    cdef int n = n0 * n1 * n2 * n3
    cdef int i

    cdef double* data = array.begin()

    cdef np.ndarray[double,ndim=1,mode='fortran'] res = np.empty(n,dtype="float64",order='F')

    for i in range(n):
        res[i] = data[i]

    return res.reshape((n0, n1, n2, n3), order='F')







