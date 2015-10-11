from bempp.core.utils.enum_types cimport HMatBlockType
from bempp.core.utils cimport shared_ptr
from bempp.core.utils cimport complex_double
from bempp.core.utils cimport Matrix

cdef extern from "bempp/hmat/hmatrix_data.hpp" namespace "hmat":
    cdef cppclass c_HMatrixData "hmat::HMatrixData"[VALUE]:
        HMatBlockType type() const


cdef extern from "bempp/hmat/hmatrix_dense_data.hpp" namespace "hmat":
    cdef cppclass c_HMatrixDenseData "hmat::HMatrixDenseData"[VALUE]:
        const Matrix[VALUE]& A() const
        int rows() const
        int cols() const
        int rank() const
        double memSizeKb() const

cdef extern from "bempp/hmat/hmatrix_low_rank_data.hpp" namespace "hmat":
    cdef cppclass c_HMatrixLowRankData "hmat::HMatrixLowRankData"[VALUE]:
        const Matrix[VALUE]& A() const
        const Matrix[VALUE]& B() const
        int rows() const
        int cols() const
        int rank() const
        double memSizeKb() const

cdef extern from "bempp/core/hmat/py_hmat_support.hpp" namespace "hmat":
    cdef shared_ptr[const c_HMatrixLowRankData[VALUE]] py_cast_to_low_rank_data[VALUE](const shared_ptr[const c_HMatrixData[VALUE]]&) 

    cdef shared_ptr[const c_HMatrixDenseData[VALUE]] py_cast_to_dense_data[VALUE](const shared_ptr[const c_HMatrixData[VALUE]]&) 


cdef class HMatrixDataBase:
    cdef _type_to_string(self, HMatBlockType block_type)

cdef class HMatrixData(HMatrixDataBase):
    cdef object _dtype
    cdef shared_ptr[const c_HMatrixData[double]] impl_float64_
    cdef shared_ptr[const c_HMatrixData[complex_double]] impl_complex128_

cdef class HMatrixLowRankData(HMatrixDataBase):
    cdef object _dtype
    cdef shared_ptr[const c_HMatrixLowRankData[double]] impl_float64_
    cdef shared_ptr[const c_HMatrixLowRankData[complex_double]] impl_complex128_

cdef class HMatrixDenseData(HMatrixDataBase):
    cdef object _dtype
    cdef shared_ptr[const c_HMatrixDenseData[double]] impl_float64_
    cdef shared_ptr[const c_HMatrixDenseData[complex_double]] impl_complex128_

cdef HMatrixLowRankData down_cast_to_low_rank_data(HMatrixData data)
cdef HMatrixDenseData down_cast_to_dense_data(HMatrixData data)




