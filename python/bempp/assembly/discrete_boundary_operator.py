"""This modules contains the data structures for assembled boundary operators."""

class DiscreteBoundaryOperator(object):
    """The base class for discretized boundary operators."""

    def __init__(self, impl):
        self._impl = impl

    @property
    def shape(self):
        """Return the shape of the discrete operator."""
        return self._impl.shape

    @property
    def dtype(self):
        """Return the type of the operator."""
        return self._impl.dtype

    def as_matrix(self):
        """Return a dense matrix representation."""

        raise NotImplementedError("Method not implemented.")

    def matvec(self, vec):
        """Multiply the operator with a numpy vector or matrix x."""
        return self._impl.matvec(vec)

    def matmat(self, mat):
        """Multiply operator with the dense numpy matrix mat."""
        return self._impl.matmat(mat)

    def __call__(self, vec):
        """Apply the operator to vec."""
        return self.matmat(vec)

    def __mul__(self, vec):

        if not isinstance(self, DiscreteBoundaryOperator):
            return vec * self

        return self.matmat(vec)

    def dot(self, vec):
        """Same as self.matvec(vec)."""

        return self.__mul__(vec)

    def __repr__(self):

        rows, cols = self.shape
        dtype = 'dtype=' + str(self.dtype)
        return "{0}x{1} {2} with {3}".format(rows, cols, self.__class__.__name__, dtype)

    def transpose(self):
        """Return the transpose of the discrete operator."""
        raise NotImplementedError("Method not implemented.")

    def adjoint(self):
        """Return the adjoint of the discrete operator."""
        raise NotImplementedError("Method not implemented.")

    def conjugate(self):
        """Return the conjugate of the discrete operator."""
        raise NotImplementedError("Method not implemented.")


class GeneralNonlocalDiscreteBoundaryOperator(DiscreteBoundaryOperator):
    """Main class for the discrete form of general discrete nonlocal operators."""

    def __init__(self, impl):

        super(GeneralNonlocalDiscreteBoundaryOperator, self).__init__(impl)

    def as_matrix(self):
        """Return a dense matrix representation."""

        return self._impl.as_matrix()

    def transpose(self):
        """Return the transpose of the discrete operator."""
        return GeneralNonlocalDiscreteBoundaryOperator(self._impl.transpose())

    def adjoint(self):
        """Return the adjoint of the discrete operator."""
        return GeneralNonlocalDiscreteBoundaryOperator(self._impl.adjoint())

    def conjugate(self):
        """Return the conjugate of the discrete operator."""
        return GeneralNonlocalDiscreteBoundaryOperator(self._impl.conjugate())


class DenseDiscreteBoundaryOperator(DiscreteBoundaryOperator):
    """Main class for the discrete form of dense discretisations of nonlocal operators."""

    def __init__(self, impl):

        super(DenseDiscreteBoundaryOperator, self).__init__(impl)

    def as_matrix(self):
        """Return a dense matrix representation."""

        return self._impl

    def matvec(self, vec):
        """Multiply the operator with a numpy vector or matrix x."""
        return self._impl.dot(vec)

    def matmat(self, mat):
        """Multiply operator with the dense numpy matrix mat."""
        return self._impl.dot(mat)

    def transpose(self):
        """Return the transpose of the discrete operator."""
        return DenseDiscreteBoundaryOperator(self._impl.T)

    def adjoint(self):
        """Return the adjoint of the discrete operator."""
        return DenseDiscreteBoundaryOperator(self._impl.T.conjugate())

    def conjugate(self):
        """Return the conjugate of the discrete operator."""
        return DenseDiscreteBoundaryOperator(self._impl.conjugate())

class SparseDiscreteBoundaryOperator(DiscreteBoundaryOperator):
    """Main class for the discrete form of sparse operators."""

    def __init__(self, impl):

        super(SparseDiscreteBoundaryOperator, self).__init__(impl)

    def as_matrix(self):
        """Return a dense matrix representation."""

        return self._impl.todense()

    def transpose(self):
        """Return the transpose of the discrete operator."""
        return SparseDiscreteBoundaryOperator(self._impl.transpose())

    def adjoint(self):
        """Return the adjoint of the discrete operator."""
        return SparseDiscreteBoundaryOperator(self._impl.transpose().conjugate())

    def conjugate(self):
        """Return the conjugate of the discrete operator."""
        return SparseDiscreteBoundaryOperator(self._impl.conjugate())

    @property
    def sparse_operator(self):
        return self._impl

class InverseSparseDiscreteBoundaryOperator(DiscreteBoundaryOperator): # pylint: disable=abstract-method
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
                                 "SparseDiscreteBoundaryOperator or of type csc_matrix.")

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

            import numpy as np
            from bempp.utils import combined_type

            if vec.ndim > 1:
                res = np.zeros((self._shape[0], vec.shape[1]), 
                               combined_type(vec.dtype, self._dtype))

                for i in range(vec.shape[1]):
                    res[:, i] = self._solve_fun(vec[:, i].ravel())
            else:
                res = self._solve_fun(vec)

            return res

        @property
        def shape(self):
            """Return the shape of the inverse operator."""
            return self._shape

        @property
        def dtype(self):
            """Return the dtype."""
            return self._dtype

    def __init__(self, operator):

        from scipy.sparse.linalg import LinearOperator
        solver = InverseSparseDiscreteBoundaryOperator._Solver(operator)

        super(InverseSparseDiscreteBoundaryOperator, self).__init__(\
                LinearOperator(solver.shape, matvec=solver.solve, rmatvec=None,
                               matmat=solver.solve, dtype=solver.dtype))




