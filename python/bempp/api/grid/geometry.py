"""Define the basic geometry class."""


class Geometry(object):
    """This class returns geometry information for an entity."""

    def __init__(self, codimension, impl):

        self._codimension = codimension
        self._impl = impl

    def integration_elements(self, local_coordinates):
        """Return the integration elements associated with the given local coordinates."""

        return self._impl.integration_elements(local_coordinates)

    def normals(self, local_coordinates):
        """Return the normal directions associated with the given element."""
        if self._codimension != 0:
            raise ValueError(
                "Method can only be called for element geometries. ")

        return self._impl.normals(local_coordinates)

    def local2global(self, local_coordinates):
        """Return the global coordinates associated with the given local coordinates."""
        return self._impl.local2global(local_coordinates)

    def jacobians_transposed(self, local_coordinates):
        """Return a list of transposed Jacobians associated with the given local coordinates."""
        return self._impl.jacobians_transposed(local_coordinates)

    def jacobian_inverses_transposed(self, local_coordinates):
        """Return a list of inverse transposed Jacobians associated with the given local coordinates."""
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
