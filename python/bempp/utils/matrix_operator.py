from bempp.assembly.discrete_boundary_operator import DiscreteBoundaryOperatorBase

class MatrixOperator(DiscreteBoundaryOperatorBase):
    def __init__(self,matrix):
        self._matrix = matrix

    def transpose(self):
        return MatrixOperator(self._matrix.transpose())

    def get_matrix(self):
        return self._matrix

    def matmat(self,x):
        return MatrixOperator(self._matrix*x)

    def matvec(self,x):
        return self._matrix.dot(x)

    def get_shape(self):
        return self._matrix.shape

    def get_dtype(self):
        return self._matrix.dtype

    def as_matrix(self):
        try:
            return self._matrix.todense()
        except AttributeError:
            return self._matrix    
    dtype = property(get_dtype)
    matrix = property(get_matrix)
    shape = property(get_shape)
    T = property(transpose)
