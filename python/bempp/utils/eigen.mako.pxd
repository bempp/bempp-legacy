<%
from data_types import dtypes,scalar_cython_type
%>


from libcpp cimport bool as cbool
cimport numpy as np
from bempp.utils cimport complex_float,complex_double

cdef extern from "bempp/common/eigen_support" namespace "Bempp":
    cdef cppclass Vector[T]:
        Vector() 
        Vector(Vector[T])
        T* data()
        T& value "operator()"(int i) 
        int rows()
        int cols()

    cdef cppclass Matrix[T]:
        Vector() 
        Vector(Vector[T])
        T* data()
        T& value "operator()"(int i, int j) 
        int rows()
        int cols()

cdef extern from "Eigen/Dense" namespace "Eigen":

% for pyvalue,cyvalue in dtypes.items():
cdef np.ndarray eigen_matrix_to_np_${pyvalue}(const Matrix[${cyvalue}]& x)
cdef np.ndarray eigen_vector_to_np_${pyvalue}(const Vector[${cyvalue}]& x)
% endfor
cdef np.ndarray eigen_matrix_to_np_int(const Matrix[int]& x)


