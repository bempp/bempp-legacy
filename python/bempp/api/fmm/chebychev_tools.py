"""Define various operations related to Chebychev Polynomials"""


class ChebychevTools(object):
    """Operations related to Chebychev Polynomials."""

    def __init__(self, order):
        """Define the object with a given order for the Chebychev polynomials."""
        from bempp.core.fmm import chebychev_tools
        self._impl = chebychev_tools.ChebychevTools(order)

    def chebychev_nodes(self):
        """Return the Chebychev nodes on the interval [-1, 1]."""

        return self._impl.chebychev_nodes()

    def evaluate_interpolation_polynomial(self, weights, evaluation_points):
        """Evaluate an interp. polynomial with given weights at the given points."""
        return self._impl.evaluate_interpolation_polynomial(
                weights, evaluation_points)

    def interpolate_to_children(self, parent_length, child_length):
        """
        Return interpolation matrix to children.

        This method returns the (2 * (order + 1) x (order + 1)) matrix, which
        maps interpolation weights at the parent level to interpolation weights
        at two children. The parent length is 'parent_length' and the child length
        is 'child_length'. The two children intervals can overlap. Hence,
        child_length > .5 * parent_length is allowed.

        """
        return self._impl.interpolate_to_children(parent_length, child_length)

    @property
    def chebychev_values(self):
        """
        Evaluate the Cheb. Polynomials at the interpolation nodes.

        Return a Numpy matrix of size (order + 1) x (order + 1), where the element
        at position (i, j) contains the jth Chebychev Polynomial evaluated
        at the ith Chebychev Point.
        """

        return self._impl.chebychev_values()
