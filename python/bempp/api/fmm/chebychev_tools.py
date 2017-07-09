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
