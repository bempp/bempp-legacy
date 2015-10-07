"""Define the basic geometry class."""

class Geometry(object):
    """This class returns geometry information for an entity."""

    def __init__(self, codimension, impl):

        self._codimension = codimension
        self._impl = impl

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



