"""This modules contains the data structures for assembled boundary operators."""

from bempp.utils.linear_operator import LinearOperator as _LinearOperator
from bempp.utils.linear_operator import MatrixLinearOperator as _MatrixLinearOperator
import numpy as _np



class GeneralNonlocalDiscreteBoundaryOperator(_LinearOperator):
    """Main class for the discrete form of general discrete nonlocal operators."""

    def __init__(self, impl):

        super(GeneralNonlocalDiscreteBoundaryOperator, self).__init__(shape=impl.shape, dtype=impl.dtype)

        self._impl = impl

    def _matvec(self, vec): # pylint: disable=method-hidden
        """Implements matrix-vector product."""

        return self._impl.matvec(vec)

    def _matmat(self, vec): # pylint: disable=method-hidden

        return self._impl.matmat(vec)

    def _adjoint(self):
        """Return the adjoint of the discrete operator."""
        return GeneralNonlocalDiscreteBoundaryOperator(self._impl.adjoint())

    def _transpose(self):
        """Return the transposed operator."""
        return GeneralNonlocalDiscreteBoundaryOperator(self._impl.transpose())

class DenseDiscreteBoundaryOperator(_MatrixLinearOperator): # pylint: disable=too-few-public-methods
    """Main class for the discrete form of dense discretisations of nonlocal operators."""

    def __init__(self, impl): # pylint: disable=super-on-old-class

        super(DenseDiscreteBoundaryOperator, self).__init__(impl)

    def __add__(self, other): # pylint: disable=super-on-old-class

        if isinstance(other, DenseDiscreteBoundaryOperator):
            return DenseDiscreteBoundaryOperator(self.A + other.A) # pylint: disable=no-member
        else:
            return super(DenseDiscreteBoundaryOperator, self).__add__(other)

    def __neg__(self):
        return DenseDiscreteBoundaryOperator(-self.A)

    def __mul__(self, other): #pylint: disable=super-on-old-class

        return self.dot(other)

    def dot(self, other):

        if isinstance(other, DenseDiscreteBoundaryOperator):
            return DenseDiscreteBoundaryOperator(self.A.dot(other.A)) #pylint: disable=no-member

        if _np.isscalar(other):
            return DenseDiscreteBoundaryOperator(self.A * other)

        return super(DenseDiscreteBoundaryOperator, self).dot(other)

    def __rmul__(self, other):

        if _np.isscalar(other):
            return DenseDiscreteBoundaryOperator(self.A * other)
        else:
            return NotImplemented

    def _transpose(self):
        """Transpose of the operator."""

        return DenseDiscreteBoundaryOperator(self.A.T)


class SparseDiscreteBoundaryOperator(_LinearOperator):
    """Main class for the discrete form of sparse operators."""

    def __init__(self, impl):

        super(SparseDiscreteBoundaryOperator, self).__init__(dtype=impl.dtype, shape=impl.shape)

        self._impl = impl

    def _matvec(self, vec):
        """Multiply the operator with a numpy vector or matrix x."""

        if self.dtype == 'float64' and _np.iscomplexobj(vec):
            return self._impl * _np.real(vec) + 1j * (self._impl * _np.imag(vec))

        return self._impl * vec

    def _matmat(self, mat):
        """Multiply operator with the dense numpy matrix mat."""
        return self._matvec(mat)

    def _transpose(self):
        """Return the transpose of the discrete operator."""
        return SparseDiscreteBoundaryOperator(self._impl.transpose())

    def _adjoint(self):
        """Return the adjoint of the discrete operator."""
        return SparseDiscreteBoundaryOperator(self._impl.transpose().conjugate())

    def __add__(self, other):

        if isinstance(other, SparseDiscreteBoundaryOperator):
            return SparseDiscreteBoundaryOperator(self.sparse_operator + other.sparse_operator)
        else:
            return super(SparseDiscreteBoundaryOperator, self).__add__(other)

    def __neg__(self):
        return SparseDiscreteBoundaryOperator(-self.sparse_operator)

    def __mul__(self, other):

        return self.dot(other)

    def dot(self, other):

        if isinstance(other, SparseDiscreteBoundaryOperator):
            return SparseDiscreteBoundaryOperator(self.sparse_operator * other.sparse_operator)

        if _np.isscalar(other):
            return SparseDiscreteBoundaryOperator(self.sparse_operator * other)

        return super(SparseDiscreteBoundaryOperator, self).dot(other)

    def __rmul__(self, other):

        if _np.isscalar(other):
            return SparseDiscreteBoundaryOperator(self.sparse_operator * other)
        else:
            return NotImplemented

    @property
    def sparse_operator(self):
        """Return the underlying Scipy sparse matrix."""
        return self._impl


class InverseSparseDiscreteBoundaryOperator(_LinearOperator):
    """Apply the (pseudo-)inverse of a sparse operator."""

    class _Solver(object): # pylint: disable=too-few-public-methods
        """Actual solver class."""

        def __init__(self, operator):

            from scipy.sparse import csc_matrix

            if isinstance(operator, SparseDiscreteBoundaryOperator):
                mat = operator.sparse_operator
            elif isinstance(operator, csc_matrix):
                mat = operator
            else:
                raise ValueError("op must be either of type " +
                                 "SparseDiscreteBoundaryOperator or of type csc_matrix. Actual type: " +
                                 str(type(operator)))

            from scipy.sparse.linalg import splu
            self._solve_fun = None
            self._shape = (mat.shape[1], mat.shape[0])
            self._dtype = mat.dtype

            if mat.shape[0] == mat.shape[1]:
                # Square matrix case
                solver = splu(mat)
                self._solve_fun = solver.solve
            elif mat.shape[0] > mat.shape[1]:
                # Thin matrix case
                mat_hermitian = mat.conjugate().transpose()
                solver = splu((mat_hermitian*mat).tocsc())
                self._solve_fun = lambda x: solver.solve(mat_hermitian*x)
            else:
                # Thick matrix case

                mat_hermitian = mat.conjugate().transpose()
                solver = splu((mat*mat_hermitian).tocsc())
                self._solve_fun = lambda x: mat_hermitian*solver.solve(x)

        def solve(self, vec):
            """Solve with right-hand side vec."""

            if self._dtype == 'float64' and _np.iscomplexobj(vec):
                return self.solve(_np.real(vec))+1j*self.solve(_np.imag(vec))

            result = self._solve_fun(vec.squeeze())

            if vec.ndim > 1:
                return result.reshape(self.shape[0], 1)
            else:
                return result

        @property
        def shape(self):
            """Return the shape of the inverse operator."""
            return self._shape

        @property
        def dtype(self):
            """Return the dtype."""
            return self._dtype

    def __init__(self, operator):

        self._solver = InverseSparseDiscreteBoundaryOperator._Solver(operator)

        super(InverseSparseDiscreteBoundaryOperator, self).__init__(\
                dtype=self._solver.dtype, shape=self._solver.shape)

    def _matvec(self, vec): #pylint: disable=method-hidden
        """Implemententation of matvec."""

        return self._solver.solve(vec)

class ZeroDiscreteBoundaryOperator(_LinearOperator):
    """A discrete operator that represents a zero operator."""

    def __init__(self, rows, columns):

        super(ZeroDiscreteBoundaryOperator, self).__init__(dtype=_np.dtype('float64'),
                                                           shape=(rows, columns))
    def _matvec(self, x):

        if x.ndim > 1:
            return _np.zeros((self.shape[0], x.shape[1]), dtype='float64')
        else:
            return _np.zeros(self.shape[0], dtype='float64')

def as_matrix(operator):
    """Return a representation of a discrete linear operator as a dense numpy matrix."""

    from numpy import eye
    cols = operator.shape[1]
    return operator * eye(cols, cols)
