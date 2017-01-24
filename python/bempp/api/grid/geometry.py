"""Define the basic geometry class."""


class Geometry(object):
    """
    This class returns geometry information for an entity.

    It is modeled onto the Dune Geometry class. A geometry
    describes the map f from the reference element onto the
    actual element and its derivative.

    Most methods take a vector of local coordinates. The reference
    triangle in BEM++ has the [x, y, z] coordinates [0, 0, 0], [1, 0, 0],
    and [0, 1, 0]. A local coordinate is given with respect to this
    reference triangle in the [x, y] plane. Hence, the vector of local
    coordinates for the three corners are

    [[0, 1, 0],
     [0, 0, 1]].

    Here, each column corresponds to one local coordinate. The
    barycenter of the triangle has the local coordinates [1./3, 1./3].
    """

    def __init__(self, codimension, impl):
        """Object constructor. Should not be called by the user."""
        self._codimension = codimension
        self._impl = impl

    def integration_elements(self, local_coordinates):
        r"""
        Return the integration elements at the given local coordinates.

        An integration element :math:`\mu` at a point :math:`p` is defined as
        :math:`\mu = \sqrt{|\det(Jf(p)^TJf(p))|}`, where :math:`Jf(p)` is the
        Jacobian of :math:`f` at :math:`p`.

        Example
        -------
        The following retunrs the integration element at the barycenter
        of an element.
        >>> p = np.array([[1./3, 1./3]]).T
        >>> mu = geom.integration_elements(p)
        """
        return self._impl.integration_elements(local_coordinates)

    def normals(self, local_coordinates):
        """Return the normal directions associated with the given element."""
        if self._codimension != 0:
            raise ValueError(
                "Method can only be called for element geometries. ")

        return self._impl.normals(local_coordinates)

    def local2global(self, local_coordinates):
        """Return global coordinates for the given local coordinates."""
        return self._impl.local2global(local_coordinates)

    def jacobians_transposed(self, local_coordinates):
        """Return transposed Jacobians for the given local coordinates."""
        return self._impl.jacobians_transposed(local_coordinates)

    def jacobian_inverses_transposed(self, local_coordinates):
        """Return the pseudo-inverses of the transposed Jacobians."""
        return self._impl.jacobian_inverses_transposed(local_coordinates)

    @property
    def corners(self):
        """Return a (3xn) array whose columns are the corners of the entity."""
        return self._impl.corners

    @property
    def corner_count(self):
        """Return the number of corners of an entity."""
        return self._impl.corner_count

    @property
    def dim(self):
        """Return the dimension of the entity."""
        return self._impl.dim

    @property
    def dim_world(self):
        """Return the dimension of the space containing the entity."""
        return self._impl.dim_world

    @property
    def volume(self):
        """Return the volume of the entity."""
        return self._impl.volume

    @property
    def codimension(self):
        """Return the codimension of the property."""
        return self._codimension
