from bempp.core.utils.enum_types cimport HMatBlockType, dense, low_rank_ab
from cython.operator cimport dereference as deref
from bempp.core.utils cimport complex_double

cdef class HMatrixDataBase:

    property block_type:

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

        def __get__(self):
            return self._get_shape()

    property rank:

        def __get__(self):
            return self._get_rank()

    property mem_size:

        def __get__(self):
            return self._get_mem_size()

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

cdef class HMatrixDenseData(HMatrixDataBase):

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
