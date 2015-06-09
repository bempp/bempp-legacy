from bempp.assembly.discrete_boundary_operator import DiscreteBoundaryOperatorBase

class MatrixOperator(DiscreteBoundaryOperatorBase):
    def __init__(self,matrix):
        self._matrix = matrix
        self._dtype = self._matrix.dtype

    def as_matrix(self):
        return self._matrix

    def matmat(self,x):
        return MatrixOperator(self._matrix*x)

    def matvec(self,x):
        return self._matrix.dot(x)

    def get_shape(self):
        return self._matrix.shape

    shape = property(get_shape)
