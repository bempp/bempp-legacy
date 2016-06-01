"""This modules contains the data structures for assembled boundary operators."""

from scipy.sparse.linalg.interface import LinearOperator as _LinearOperator
import numpy as _np


def _call_super(obj, dtype, shape):
    """Call the correct super constructor depending on scipy version."""


class DiscreteBoundaryOperator(_LinearOperator):
    """Base class for discrete boundary operators."""

    def __new__(cls, *args, **kwargs):

        # Overwriting new because LinearOperator calls __init__
        # unnecessarily in its __new__ method causing doubly
        # called constructors (to be fixed in 0.18)

        return object.__new__(cls)

    def __init__(self, dtype, shape):

        import scipy
        if scipy.__version__ < '0.16.0':
            super(DiscreteBoundaryOperator, self).__init__(shape, self._matvec, rmatvec=self._rmatvec,
                                                           matmat=self._matmat, dtype=dtype)
        else:
            super(DiscreteBoundaryOperator, self).__init__(dtype, shape)

    def __add__(self, other):

        if isinstance(other, DiscreteBoundaryOperator):
            return DiscreteBoundaryOperatorSum(self, other)
        else:
            return super(DiscreteBoundaryOperator, self).__add__(other)

    def __mul__(self, other):

        return self.dot(other)

    def dot(self, other):

        if isinstance(other, DiscreteBoundaryOperator):
            return DiscreteBoundaryOperatorProduct(self, other)
        elif isinstance(other, _LinearOperator):
            return super(DiscreteBoundaryOperator, self).dot(other)
        elif _np.isscalar(other):
            return ScaledDiscreteBoundaryOperator(self, other)
        else:
            x = _np.asarray(other)
            if x.ndim == 1 or (x.ndim == 2 and x.shape[1] == 1):
                return self._matvec(x)
            elif x.ndim == 2:
                return self._matmat(x)
            else:
                raise ValueError("Expect a 1d or 2d array or matrix.")

    def __rmul__(self, other):

        if _np.isscalar(other):
            return self * other
        else:
            raise ValueError(
                "Cannot multiply operand of type {0} from the left.".format(type(other)))

    def __call__(self, other):

        return self.dot(other)

    def __matmul__(self, other):

        if np.isscalar(other):
            raise ValueError("Scalar operands not allowed. Use '*' instead.")

        return self.dot(other)

    def __neg__(self):

        return -1 * self

    def __sub__(self, other):

        return self.__add__(-other)


class DiscreteBoundaryOperatorSum(DiscreteBoundaryOperator):

    def __init__(self, op1, op2):

        if not isinstance(op1, DiscreteBoundaryOperator) or \
                not isinstance(op2, DiscreteBoundaryOperator):
            raise ValueError(
                "Both operators must be discrete boundary operators.")

        if op1.shape != op2.shape:
            raise ValueError(
                "Shape mismatch: {0} != {1}.".format(op1.shape, op2.shape))

        self._op1 = op1
        self._op2 = op2

        super(DiscreteBoundaryOperatorSum, self).__init__(
            _np.find_common_type([op1.dtype, op2.dtype], []),
            op1.shape)

    def _matvec(self, x):

        return self._op1.matvec(x) + self._op2.matvec(x)

    def _matmat(self, x):

        return self._op1.matmat(x) + self._op2.matmat(x)

    def _rmatvec(self, x):

        return self._op1.rmatvec(x) + self._op2.rmatvec(x)

    def _adjoint(self, x):

        return self._op1.adjoint() + self._op2.adjoint()

    def _transpose(self, x):

        return self._op1.transpose() + self._op2.transpose()


class DiscreteBoundaryOperatorProduct(DiscreteBoundaryOperator):

    def __init__(self, op1, op2):

        if not isinstance(op1, DiscreteBoundaryOperator) or \
                not isinstance(op2, DiscreteBoundaryOperator):
            raise ValueError(
                "Both operators must be discrete boundary operators.")

        if op1.shape[1] != op2.shape[0]:
            raise ValueError("Shapes {0} and {1} not compatible for matrix product.".format(
                op1.shape, op2.shape))

        self._op1 = op1
        self._op2 = op2

        super(DiscreteBoundaryOperatorProduct, self).__init__(
            _np.find_common_type([op1.dtype, op2.dtype], []),
            (op1.shape[0], op2.shape[1]))

    def _matvec(self, x):

        return self._op1.matvec(self._op2.matvec(x))

    def _matmat(self, x):

        return self._op1.matmat(self._op2.matmat(x))

    def _rmatvec(self, x):

        return self._op2.rmatvec(self._op1.rmatvec(x))

    def _adjoint(self, x):

        return self._op2.adjoint() * self._op1.adjoint()

    def _transpose(self, x):

        return self._op2.transpose() + self._op1.transpose()


class ScaledDiscreteBoundaryOperator(DiscreteBoundaryOperator):

    def __init__(self, op, alpha):

        if not isinstance(op, DiscreteBoundaryOperator):
            raise ValueError(
                "Both operators must be discrete boundary operators.")

        self._op = op
        self._alpha = alpha

        super(ScaledDiscreteBoundaryOperator, self).__init__(
            _np.find_common_type([op.dtype, _np.array([alpha]).dtype], []),
            op.shape)

    def _matvec(self, x):

        return self._alpha * self._op.matvec(x)

    def _matmat(self, x):

        return self._alpha * self._op.matmat(x)

    def _rmatvec(self, x):

        return self._alpha * self._op.rmatvec(x)

    def _adjoint(self, x):

        return self._alpha * self._op.adjoint()

    def _transpose(self, x):

        return self._alpha * self._op.transpose()


class GeneralNonlocalDiscreteBoundaryOperator(DiscreteBoundaryOperator):
    """Main class for the discrete form of general discrete nonlocal operators.

    This class derives from :class:`scipy.sparse.linalg.interface.LinearOperator`
    and thereby implements the SciPy LinearOperator protocol.

    """

    def __init__(self, impl):

        super(GeneralNonlocalDiscreteBoundaryOperator,
              self).__init__(impl.dtype, impl.shape)

        self._impl = impl

    def _matvec(self, vec):  # pylint: disable=method-hidden
        """Implements matrix-vector product."""

        return self._impl.matvec(vec)

    def _matmat(self, vec):  # pylint: disable=method-hidden

        return self._impl.matmat(vec)

    def _adjoint(self):
        """Return the adjoint of the discrete operator."""
        return GeneralNonlocalDiscreteBoundaryOperator(self._impl.adjoint())

    def _transpose(self):
        """Return the transposed operator."""
        return GeneralNonlocalDiscreteBoundaryOperator(self._impl.transpose())


class DenseDiscreteBoundaryOperator(DiscreteBoundaryOperator):  # pylint: disable=too-few-public-methods
    """Main class for the discrete form of dense discretisations of nonlocal operators.

    This class derives from :class:`scipy.sparse.linalg.interface.LinearOperator`
    and thereby implements the SciPy LinearOperator protocol.

    """

    def __init__(self, impl):  # pylint: disable=super-on-old-class

        self._impl = impl
        super(DenseDiscreteBoundaryOperator, self).__init__(
            impl.dtype, impl.shape)

    def _matvec(self, x):

        return self._matmat(x)

    def _matmat(self, x):

        return self.A.dot(x)

    def _rmatvec(self, x):

        return x.dot(self.A)

    def __add__(self, other):  # pylint: disable=super-on-old-class

        if isinstance(other, DenseDiscreteBoundaryOperator):
            return DenseDiscreteBoundaryOperator(self.A + other.A)  # pylint: disable=no-member
        else:
            return super(DenseDiscreteBoundaryOperator, self).__add__(other)

    def __neg__(self):
        return DenseDiscreteBoundaryOperator(-self.A)

    def __mul__(self, other):  # pylint: disable=super-on-old-class

        return self.dot(other)

    def dot(self, other):

        if isinstance(other, DenseDiscreteBoundaryOperator):
            return DenseDiscreteBoundaryOperator(self.A.dot(other.A))  # pylint: disable=no-member

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

    def _adjoint(self):
        """Adjoint of the operator."""

        return DenseDiscreteBoundaryOperator(self.A.conjugate().transpose())

    @property
    def A(self):
        """Return the underlying array."""

        return self._impl


class SparseDiscreteBoundaryOperator(DiscreteBoundaryOperator):
    """Main class for the discrete form of sparse operators.

    This class derives from :class:`scipy.sparse.linalg.interface.LinearOperator`
    and thereby implements the SciPy LinearOperator protocol.

    """

    def __init__(self, impl):

        super(SparseDiscreteBoundaryOperator, self).__init__(
            impl.dtype, impl.shape)

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


class InverseSparseDiscreteBoundaryOperator(DiscreteBoundaryOperator):
    """Apply the (pseudo-)inverse of a sparse operator.

    This class uses a Sparse LU-Decomposition (in the case of a square matrix)
    or a sparse normal equation to provide the application of an inverse to
    a sparse operator.

    This class derives from :class:`scipy.sparse.linalg.interface.LinearOperator`
    and thereby implements the SciPy LinearOperator protocol.

    Parameters
    ----------
    operator : bempp.api.SparseDiscreteBoundaryOperator or scipy.sparse.csc_matrix
        Sparse operator to be inverted.

    """

    class _Solver(object):  # pylint: disable=too-few-public-methods
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

            from scipy.sparse.linalg import factorized
            from scipy.sparse.linalg import splu
            self._solve_fun = None
            self._shape = (mat.shape[1], mat.shape[0])
            self._dtype = mat.dtype

            import time
            import bempp.api

            bempp.api.LOGGER.info("Start computing LU (pseudo)-inverse of ({0}, {1}) matrix.".format(
                mat.shape[0], mat.shape[1]))

            start_time = time.time()
            if mat.shape[0] == mat.shape[1]:
                # Square matrix case
                self._solve_fun = factorized(mat)
            elif mat.shape[0] > mat.shape[1]:
                # Thin matrix case
                mat_hermitian = mat.conjugate().transpose()
                solver = splu((mat_hermitian * mat).tocsc())
                self._solve_fun = lambda x: solver.solve(mat_hermitian * x)
            else:
                # Thick matrix case

                mat_hermitian = mat.conjugate().transpose()
                solver = splu((mat * mat_hermitian).tocsc())
                self._solve_fun = lambda x: mat_hermitian * solver.solve(x)

            end_time = time.time()
            bempp.api.LOGGER.info("Finished computation of inverse in {0:.2E} seconds.".format(
                end_time - start_time))

        def solve(self, vec):
            """Solve with right-hand side vec."""

            if self._dtype == 'float64' and _np.iscomplexobj(vec):
                return self.solve(_np.real(vec)) + 1j * self.solve(_np.imag(vec))

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
        super(DiscreteBoundaryOperator, self).__init__(
            self._solver.dtype, self._solver.shape)

    def _matvec(self, vec):  # pylint: disable=method-hidden
        """Implemententation of matvec."""

        return self._solver.solve(vec)


class ZeroDiscreteBoundaryOperator(DiscreteBoundaryOperator):
    """A discrete operator that represents a zero operator.

    This class derives from :class:`scipy.sparse.linalg.interface.LinearOperator`
    and thereby implements the SciPy LinearOperator protocol.

    Parameters
    ----------
    rows : int
        The number of rows in the operator.
    columns : int
        The number of columns in the operator.

    """

    def __init__(self, rows, columns):

        super(ZeroDiscreteBoundaryOperator, self).__init__(_np.dtype('float64'),
                                                           (rows, columns))

    def _matvec(self, x):

        if x.ndim > 1:
            return _np.zeros((self.shape[0], x.shape[1]), dtype='float64')
        else:
            return _np.zeros(self.shape[0], dtype='float64')


class DiscreteRankOneOperator(DiscreteBoundaryOperator):
    """Creates a discrete rank one operator.

    This class represents a rank one operator given
    by column * row, where column is column is
    interpreted as a (m, 1) array and row as
    a (1, n) array.

    Parameters
    ----------
    column : np.array
        A column vector
    row : np.array
        A row vector

    """

    def __init__(self, column, row):

        if row.dtype == 'complex128' or column.dtype == 'complex128':
            dtype = 'complex128'
        else:
            dtype = 'float64'

        self._row = row.ravel()
        self._column = column.ravel()

        shape = (len(self._column), len(self._row))

        super(DiscreteRankOneOperator, self).__init__(dtype,
                                                      shape)

    def _matvec(self, x):

        import numpy as np

        if x.ndim > 1:
            return np.outer(self._column, np.dot(self._row, x))
        else:
            return self._column * np.dot(self._row, x)

    def _transpose(self, x):

        return DiscreteRankOneOperator(row, column)

    def _adjoint(self, x):

        return DiscreteRankOneOperator(row.conjugate(), column.conjugate())


def as_matrix(operator):
    """Return a representation of a discrete linear operator as a dense numpy matrix.

    Parameters
    ----------
    operator : scipy.sparse.linalg.interface.LinearOperator
        The linear operator to be converted into a dense matrix.


    Notes
    -----
    Note that this function may be slow depending on how the original discrete operator was stored.
    In the case of a dense assembly simple the underlying NumPy matrix is returned. Otherwise,
    the operator needs to be converted to an array, which can take a long time.

    """

    from numpy import eye
    cols = operator.shape[1]
    if isinstance(operator, DenseDiscreteBoundaryOperator):
        return operator.A
    elif isinstance(operator, SparseDiscreteBoundaryOperator):
        return operator.sparse_operator
    else:
        return operator * eye(cols, cols)
