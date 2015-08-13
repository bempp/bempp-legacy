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
        return self.matvec(mat)

    def __call__(self, vec):
        """Apply the operator to vec."""
        return self.matvec(vec)

    def __mul__(self, vec):

        if not isinstance(self, DiscreteBoundaryOperator):
            return vec * self

        return self.matvec(vec)

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


