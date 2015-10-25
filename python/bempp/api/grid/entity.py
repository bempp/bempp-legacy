"""Define the basic entity class."""
from bempp.api.grid.geometry import Geometry as _Geometry

class Entity(object):
    """This class provides an interface to grid entities."""

    def __init__(self, codimension, impl):

        self._codimension = codimension
        self._impl = impl


    def sub_entity_iterator(self, codimension):

        if self.codimension != 0:
            raise Exception("Error. Codimension is {0}, but method can only be called for entities of codimension 0.".format(self.codimension))

        from .entity_iterator import EntityIterator
        return EntityIterator(codimension,
                              self._impl.sub_entity_iterator(codimension))


    @property
    def domain(self):
        """Return the domain index of the entity."""
        if self.codimension != 0:
            raise Exception("Error. Codimension is {0}, but attribute only exists for entities of codimension 0.".format(self.codimension))

        return self._impl.domain

    @property
    def is_leaf(self):
        """Return true if the entity has no descendents."""
        if self.codimension != 0:
            raise Exception("Error. Codimension is {0}, but attribute only exists for entities of codimension 0.".format(self.codimension))

        return self._impl.is_leaf

    @property
    def has_father(self):
        """Return true if the entity has a father entity."""
        if self.codimension != 0:
            raise Exception("Error. Codimension is {0}, but attribute only exists for entities of codimension 0.".format(self.codimension))

        return self._impl.has_father

    @property
    def father(self):
        """Return father entity."""
        if self.codimension != 0:
            raise Exception("Error. Codimension is {0}, but attribute only exists for entities of codimension 0.".format(self.codimension))

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



