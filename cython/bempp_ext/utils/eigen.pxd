from libcpp cimport bool as cbool
from bempp_ext.utils.complex cimport complex_double
cimport numpy as np

cdef extern from "bempp/common/eigen_support.hpp" namespace "Bempp":
    cdef cppclass Vector[T]:
        Vector()
        Vector(int)
        Vector(Vector[T])
        void resize(int)
        T* data()
        T& value "operator()"(int i) 
        int rows()
        int cols()

    cdef cppclass Matrix[T]:
        Matrix()
        Matrix(int,int)
        Matrix(Matrix[T])
        void resize(int,int)
        T* data()
        T& value "operator()"(int i, int j) 
        int rows()
        int cols()

cdef extern from "bempp_ext/utils/py_utils.hpp" namespace "Bempp":
    cdef Vector[T] copy_buf_to_vec[T](T* buf, int n)
    cdef Matrix[T] copy_buf_to_mat[T](T* buf, int m, int n)

cdef np.ndarray eigen_matrix_to_np_float64(const Matrix[double]& x)
cdef np.ndarray eigen_matrix_to_np_complex128(const Matrix[complex_double]& x)

cdef np.ndarray eigen_vector_to_np_float64(const Vector[double]& x)
cdef np.ndarray eigen_vector_to_np_complex128(const Vector[complex_double]& x)

cdef Matrix[double] np_to_eigen_matrix_float64(np.ndarray x)
cdef Matrix[complex_double] np_to_eigen_matrix_complex128(np.ndarray x)

cdef Vector[double] np_to_eigen_vector_float64(np.ndarray x)
cdef Vector[complex_double] np_to_eigen_vector_complex128(np.ndarray x)

cdef np.ndarray eigen_matrix_to_np_int(const Matrix[int]& x)
cdef np.ndarray eigen_vector_to_np_int(const Vector[int]& x)




