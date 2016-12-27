"""
Definition of shapesets.

A shapeset is the set of basis functions on a given reference element
associated with an actual element.

"""

class Shapeset(object):
    """Access the data of the basis functions on an element."""

    def __init__(self, impl):
        self._impl = impl

    def evaluate(self, points, dof, values=True, derivatives=False):
        """
        Evaluate the shapeset.

        Parameters
        ----------
        points : numpy.ndarray
            A 3 x n array of evaluation points.
        dof : int
            The local index of the dof on the shapeset. To evaluate on
            all dofs use bempp.api.ALL.

        Returns
        -------
        An ndarray object. The last dimension corresponds to the number of
        points, the second last dimension to the number of dofs, and the
        first one or two dimensions to the function value, respectively
        derivative matrix. If values and derivatives are both True the
        return value is a tuple consisting of values and derivatives.

        """
        return self._impl.evaluate(points, dof, values, derivatives)


    @property
    def order(self):
        """Order of the shapeset."""
        return self._impl.order

    @property
    def size(self):
        """Number of basis functions in the shapeset."""
        return self._impl.size
