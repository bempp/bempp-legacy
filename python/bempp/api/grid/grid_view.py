"""Implementation of a view onto a grid."""


class GridView(object):
    """Provides interfaces to query the entities of a grid."""

    def __init__(self, impl):
        self._impl = impl

    def entity_count(self, codimension):
        """Return the number of entities of a given codimension."""
        return self._impl.entity_count(codimension)

    def index_set(self):
        """Return an IndexSet object for the GridView."""
        from .index_set import IndexSet
        return IndexSet(self._impl.index_set())

    def entity_iterator(self, codimension):
        """Return an entity iterator for a given codimension."""
        from .entity_iterator import EntityIterator
        return EntityIterator(codimension,
                              self._impl.entity_iterator(codimension))

    @property
    def dim(self):
        """Return the dimension of the grid."""
        return self._impl.dim

    @property
    def dim_world(self):
        """Return the dimension of the space containing the grid."""
        return self._impl.dim_world

    @property
    def vertices(self):
        """Return a (3 x n_vertices) array with the vertices of the grid."""
        return self._impl.vertices

    @property
    def elements(self):
        """Return a (3 x n_elements) array with the elements of the grid."""
        return self._impl.elements

    @property
    def domain_indices(self):
        """Return a list of domain indices."""
        return self._impl.domain_indices
