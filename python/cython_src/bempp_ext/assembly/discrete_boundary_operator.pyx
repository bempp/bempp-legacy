from bempp_ext.utils cimport eigen_matrix_to_np_float64
from bempp_ext.utils cimport eigen_matrix_to_np_complex128
from bempp_ext.utils cimport enum_types as enums
from libcpp cimport bool as cbool
from cython.operator cimport dereference as deref

import numpy as np
cimport numpy as np


cdef class RealDiscreteBoundaryOperator:

    def __cinit__(self):
        self.transpose_mode = enums.no_transpose

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()

    property dtype:
        def __get__(self):
            return np.dtype('float64')

    property transpose_mode:
        def __get__(self):
            if self.transpose_mode == enums.no_transpose:
                return 'no_transpose'
            if self.transpose_mode == enums.conjugate:
                return 'conjugate'
            if self.transpose_mode == enums.transpose:
                return 'transpose'
            if self.transpose_mode == enums.conjugate_transpose:
                return 'conjugate_transpose'

        def __set__(self, value):
            self.transpose_mode = enums.transposition_mode(value)

    property shape:
        def __get__(self):
            cdef unsigned int rows
            cdef unsigned int cols

            rows = deref(self.impl_).rowCount()
            cols = deref(self.impl_).columnCount()

            if (self.transpose_mode==enums.no_transpose or
                    self.transpose_mode==enums.conjugate):
                return (rows,cols)
            else:
                return (cols,rows)

    def as_matrix(self):
        cdef Matrix[double] mat_data = deref(self.impl_).asMatrix()
        if (self.transpose_mode==enums.no_transpose):
            return eigen_matrix_to_np_float64(mat_data)
        if (self.transpose_mode==enums.conjugate):
            return np.conjugate(eigen_matrix_to_np_float64(mat_data))
        if (self.transpose_mode==enums.transpose):
            return np.transpose(eigen_matrix_to_np_float64(mat_data))
        if (self.transpose_mode==enums.conjugate_transpose):
            return np.conjugate(
                    np.transpose(eigen_matrix_to_np_float64(mat_data)))

    def matvec(self,np.ndarray x):

        if np.iscomplexobj(x):
            return self*np.real(x)+1j*(self*np.imag(x))

        cdef np.ndarray x_in
        cdef np.ndarray y_inout

        cdef int rows = self.shape[0]
        cdef int cols = self.shape[1]

        cdef cbool is_reshaped = False

        if (x.ndim==1):
            x_in = x.reshape((-1,1),order='F').astype(self.dtype,
                    order='F',casting='safe',copy=False)
            is_reshaped = True
        elif (x.ndim==2):
            x_in = x.astype(self.dtype,order='F',casting='safe',copy=False)
        else:
            raise ValueError('x must have at most two dimensions')

        if (x_in.shape[0]!=self.shape[1]):
            raise ValueError("Wrong dimensions.")

        y = deref(self.impl_).apply(self.transpose_mode,x_in)

        if is_reshaped:
            y = y.ravel()

        return y

    def matmat(self, np.ndarray x):
        return self.matvec(x)

    def transpose(self):

        cdef RealDiscreteBoundaryOperator op = RealDiscreteBoundaryOperator()
        op.impl_ = self.impl_
        op.transpose_mode = enums.compute_transpose_mode(self.transpose_mode,
                                                          enums.transpose)
        return op

    def adjoint(self):

        cdef RealDiscreteBoundaryOperator op = RealDiscreteBoundaryOperator()
        op.impl_ = self.impl_
        op.transpose_mode = enums.compute_transpose_mode(self.transpose_mode,
                                                          enums.adjoint)
        return op

    def conjugate(self):

        cdef RealDiscreteBoundaryOperator op = RealDiscreteBoundaryOperator()
        op.impl_ = self.impl_
        op.transpose_mode = enums.compute_transpose_mode(self.transpose_mode,
                                                          enums.conjugate)
        return op

cdef class ComplexDiscreteBoundaryOperator:

    def __cinit__(self):
        self.transpose_mode = enums.no_transpose

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()

    property dtype:
        def __get__(self):
            return np.dtype('complex128')

    property transpose_mode:
        def __get__(self):
            if self.transpose_mode == enums.no_transpose:
                return 'no_transpose'
            if self.transpose_mode == enums.conjugate:
                return 'conjugate'
            if self.transpose_mode == enums.transpose:
                return 'transpose'
            if self.transpose_mode == enums.conjugate_transpose:
                return 'conjugate_transpose'

        def __set__(self, value):
            self.transpose_mode = enums.transposition_mode(value)

    property shape:
        def __get__(self):
            cdef unsigned int rows
            cdef unsigned int cols

            rows = deref(self.impl_).rowCount()
            cols = deref(self.impl_).columnCount()

            if (self.transpose_mode==enums.no_transpose or
                    self.transpose_mode==enums.conjugate):
                return (rows,cols)
            else:
                return (cols,rows)

    def as_matrix(self):
        cdef Matrix[complex_double] mat_data = deref(self.impl_).asMatrix()
        if (self.transpose_mode==enums.no_transpose):
            return eigen_matrix_to_np_complex128(mat_data)
        if (self.transpose_mode==enums.conjugate):
            return np.conjugate(eigen_matrix_to_np_complex128(mat_data))
        if (self.transpose_mode==enums.transpose):
            return np.transpose(eigen_matrix_to_np_complex128(mat_data))
        if (self.transpose_mode==enums.conjugate_transpose):
            return np.conjugate(
                    np.transpose(eigen_matrix_to_np_complex128(mat_data)))

    def matvec(self,np.ndarray x):

        cdef np.ndarray x_in
        cdef np.ndarray y_inout

        cdef int rows = self.shape[0]
        cdef int cols = self.shape[1]

        cdef cbool is_reshaped = False

        if (x.ndim==1):
            x_in = x.reshape((-1,1),order='F').astype(self.dtype,
                    order='F',casting='safe',copy=False)
            is_reshaped = True
        elif (x.ndim==2):
            x_in = x.astype(self.dtype,order='F',casting='safe',copy=False)
        else:
            raise ValueError('x must have at most two dimensions')

        if (x_in.shape[0]!=self.shape[1]):
            raise ValueError("Wrong dimensions.")

        y = deref(self.impl_).apply(self.transpose_mode,x_in)

        if is_reshaped:
            y = y.ravel()

        return y

    def matmat(self, np.ndarray x):
        return self.matvec(x)

    def transpose(self):

        cdef ComplexDiscreteBoundaryOperator op = ComplexDiscreteBoundaryOperator()
        op.impl_ = self.impl_
        op.transpose_mode = enums.compute_transpose_mode(self.transpose_mode,
                                                          enums.transpose)
        return op

    def adjoint(self):

        cdef ComplexDiscreteBoundaryOperator op = ComplexDiscreteBoundaryOperator()
        op.impl_ = self.impl_
        op.transpose_mode = enums.compute_transpose_mode(self.transpose_mode,
                                                          enums.adjoint)
        return op

    def conjugate(self):

        cdef ComplexDiscreteBoundaryOperator op = ComplexDiscreteBoundaryOperator()
        op.impl_ = self.impl_
        op.transpose_mode = enums.compute_transpose_mode(self.transpose_mode,
                                                          enums.conjugate)
        return op

def convert_to_sparse(RealDiscreteBoundaryOperator op):
    """Convert sparse discrete boundary operator to CSC matrix."""
    return c_convert_to_sparse(op.impl_)

def convert_to_dense(op):
    """Convert dense discrete boundary operator to numpy array."""
    
    if isinstance(op, RealDiscreteBoundaryOperator):
        return c_convert_to_dense_double((<RealDiscreteBoundaryOperator>op).impl_)
    elif isinstance(op, ComplexDiscreteBoundaryOperator):
        return c_convert_to_dense_complex((<ComplexDiscreteBoundaryOperator>op).impl_)
    else:
        raise ValueError("Type not supported.")



