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

    def child_interpolation_matrix(self, ratio):
        """
        Return interpolation matrix to children.

        This method returns the (2 * (order + 1) x (order + 1)) matrix, which
        maps interpolation weights at the parent level to interpolation weights
        at two children. 'ratio' is the ratio of the width of one child to the
        width of the parent. The children can overlap. Hence, the ratio can be
        larger than .5

        """
        return self._impl.child_interpolation_matrix(ratio)

    def derivative_weights(self, weights):
        """Interpolate the derivative on the interval [-1, 1]."""
        return self._impl.derivative_weights(weights)

    def derivative_weights_3d(self, weights, direction ):
        """Interpolate the derivative on the cube [-1, 1]^3."""
        return self._impl.derivative_weights_3d(weights, direction)

