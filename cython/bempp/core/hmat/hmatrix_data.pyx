# cython: embedsignature=True

from bempp.core.utils.enum_types cimport HMatBlockType, dense, low_rank_ab
from cython.operator cimport dereference as deref
from bempp.core.utils cimport complex_double
from bempp.core.utils cimport eigen_matrix_to_np_float64
from bempp.core.utils cimport eigen_matrix_to_np_complex128

cdef class HMatrixDataBase:
    """Base class for H-Matrix data."""

    property block_type:
        """Returns the type of the data block (`dense` or `low_rank_ab`). """

        def __get__(self):
            return self._get_type()

    def _get_type(self):
        pass

    cdef _type_to_string(self, HMatBlockType block_type):

        if block_type == dense:
            return 'dense'
        if block_type == low_rank_ab:
            return 'low_rank_ab'

        raise ValueError("Unsupported block type.")

    property shape:
        """Returns the shape of the data block."""

        def __get__(self):
            return self._get_shape()

    property rank:
        """Returns the rank of the data."""

        def __get__(self):
            return self._get_rank()

    property mem_size:
        """Returns the memory size in kb of the data."""

        def __get__(self):
            return self._get_mem_size()

    property dtype:
        """Returns the data type."""

        def __get__(self):
            return self._dtype

cdef class HMatrixData(HMatrixDataBase):

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):

        self.impl_float64_.reset()
        self.impl_complex128_.reset()

    def _get_type(self):

        if self._dtype=='float64':
            return self._type_to_string(deref(self.impl_float64_).type())

        if self._dtype=='complex128':
            return self._type_to_string(deref(self.impl_complex128_).type())


        
cdef class HMatrixLowRankData(HMatrixDataBase):
    """Interface to low-rank H-Matrix data.

    This class provides an interface to low-rank H-Matrix data.
    Algebraically, the stored matrix D has the form D = A * B,
    where A is an (m, k) matrix and B a (k, n) matrix. Here,
    k is the rank, m is the number of rows and n is the number
    of columns.

    """




    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        pass

    def _get_type(self):
        return 'low_rank_ab'

    def _get_shape(self):
        cdef int rows
        cdef int cols

        if self._dtype == 'float64':
            rows = deref(self.impl_float64_).rows()
            cols = deref(self.impl_float64_).cols()
            return (rows,cols)

        if self._dtype == 'complex128':
            rows = deref(self.impl_complex128_).rows()
            cols = deref(self.impl_complex128_).cols()
            return (rows,cols)

    def _get_rank(self):

        if self._dtype == 'float64':
            return deref(self.impl_float64_).rank()

        if self._dtype == 'complex128':
            return deref(self.impl_complex128_).rank()
        raise ValueError("Unsupported dtype.")

    def _get_mem_size(self):

        if self._dtype == 'float64':
            return deref(self.impl_float64_).memSizeKb()

        if self._dtype == 'complex128':
            return deref(self.impl_complex128_).memSizeKb()
        raise ValueError("Unsupported dtype.")

    property A:
        """Returns the matrix A as NumPy array."""

        def __get__(self):

            if self.dtype == 'float64':
                return eigen_matrix_to_np_float64(deref(self.impl_float64_).A())
            else:
                return eigen_matrix_to_np_complex128(deref(self.impl_complex128_).A())
    property B:
        """Returns the matrix B as NumPy array."""

        def __get__(self):

            if self.dtype == 'float64':
                return eigen_matrix_to_np_float64(deref(self.impl_float64_).B())
            else:
                return eigen_matrix_to_np_complex128(deref(self.impl_complex128_).B())


cdef class HMatrixDenseData(HMatrixDataBase):
    """Interface to dense H-Matrix data.

    This class provides an interface to dense H-Matrix data, that
    is internally the data is stored as a dense matrix A.

    """



    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        pass

    def _get_type(self):
        return 'dense'

    def _get_shape(self):
        cdef int rows
        cdef int cols

        if self._dtype == 'float64':
            rows = deref(self.impl_float64_).rows()
            cols = deref(self.impl_float64_).cols()
            return (rows,cols)

        if self._dtype == 'complex128':
            rows = deref(self.impl_complex128_).rows()
            cols = deref(self.impl_complex128_).cols()
            return (rows,cols)

        raise ValueError("Unsupported dtype.")

    def _get_rank(self):

        if self._dtype == 'float64':
            return deref(self.impl_float64_).rank()

        if self._dtype == 'complex128':
            return deref(self.impl_complex128_).rank()
        raise ValueError("Unsupported dtype.")

    def _get_mem_size(self):

        if self._dtype == 'float64':
            return deref(self.impl_float64_).memSizeKb()

        if self._dtype == 'complex128':
            return deref(self.impl_complex128_).memSizeKb()
        raise ValueError("Unsupported dtype.")

    property A:
        """Return the dense matrix."""

        def __get__(self):

            if self.dtype == 'float64':
                return eigen_matrix_to_np_float64(deref(self.impl_float64_).A())
            else:
                return eigen_matrix_to_np_complex128(deref(self.impl_complex128_).A())

cdef HMatrixLowRankData down_cast_to_low_rank_data(HMatrixData data):

    cdef HMatrixLowRankData result = HMatrixLowRankData()

    if data._dtype == 'float64':
        result._dtype = data._dtype
        result.impl_float64_.assign(
                py_cast_to_low_rank_data[double](data.impl_float64_))
    elif data._dtype == 'complex128':
        result._dtype = data._dtype
        result.impl_complex128_.assign(
                py_cast_to_low_rank_data[complex_double](data.impl_complex128_))
    else:
        raise ValueError("Unsupported dtype.")

    return result

cdef HMatrixDenseData down_cast_to_dense_data(HMatrixData data):

    cdef HMatrixDenseData result = HMatrixDenseData()

    if data._dtype == 'float64':
        result._dtype = data._dtype
        result.impl_float64_.assign(
                py_cast_to_dense_data[double](data.impl_float64_))
    elif data._dtype == 'complex128':
        result._dtype = data._dtype
        result.impl_complex128_.assign(
                py_cast_to_dense_data[complex_double](data.impl_complex128_))
    else:
        raise ValueError("Unsupported dtype.")

    return result
