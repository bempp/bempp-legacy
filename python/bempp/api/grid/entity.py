"""Define the basic entity class."""
from bempp.api.grid.geometry import Geometry as _Geometry


class Entity(object):
    """
    This class provides an interface to grid entities.

    Entities in BEM++ are vertices, edges or elements. The type of an
    entity depends on its codimension. Elements have codimension 0,
    edges have codimension 1, and vertices have codimension 2 with
    respect to the grid dimension.

    """

    def __init__(self, codimension, impl):
        """Constructor. Should not be called by the user."""
        self._codimension = codimension
        self._impl = impl

    def sub_entity_iterator(self, codimension):
        """
        Return an iterator for subentities of given codimension.

        Note that this method can only be called for entities of
        codimension 0.

        Example
        -------
        The following iterates through the vertices of an element and
        prints their coordinates.

        >>> for vertex in element.sub_entity_iterator(2):
        >>>     print(vertex.geometry.corners)

        """
        import bempp.api.grid.entity_iterator as _entity_iterator
        if self.codimension != 0:
            raise Exception(
                "Error. Codimension is {0}".format(self.codimension) +
                ", but method can only be called for entities of codim 0.")

        return _entity_iterator.EntityIterator(
            codimension,
            self._impl.sub_entity_iterator(codimension))

    @property
    def domain(self):
        """Return the domain index of the entity."""
        if self.codimension != 0:
            raise Exception(
                "Error. Codimension is {0}".format(self.codimension) +
                ", but attribute only exists for entities of codim 0.")

        return self._impl.domain

    @property
    def is_leaf(self):
        """Return true if the entity has no descendents."""
        if self.codimension != 0:
            raise Exception(
                "Error. Codimension is {0}".format(self.codimension) +
                ", but attribute only exists for entities of codim 0.")

        return self._impl.is_leaf

    @property
    def has_father(self):
        """Return true if the entity has a father entity."""
        if self.codimension != 0:
            raise Exception(
                "Error. Codimension is {0}".format(self.codimension) +
                ", but attribute only exists for entities of codim 0.")

        return self._impl.has_father

    @property
    def father(self):
        """Return father entity."""
        if self.codimension != 0:
            raise Exception(
                "Error. Codimension is {0}".format(self.codimension) +
                ", but attribute only exists for entities of codim 0.")

        return self._impl.father

    @property
    def level(self):
        """Return the level of the entity."""
        return self._impl.level()

    @property
    def geometry(self):
        """Return the associated geometry."""
        return _Geometry(self.codimension, self._impl.geometry)

    @property
    def codimension(self):
        """Return the codimension of the entity."""
        return self._codimension
