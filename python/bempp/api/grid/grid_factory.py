"""Interfaces to a grid factory."""


class GridFactory(object):
    """
    Implement a grid factory.

    A grid factory can be used to create arbitrary
    grids from a set of vertices and elements.

    Examples
    --------
    The following gives an example of how to grid an grid
    containing a single element using a grid factory.

    >>> factory = bempp.api.grid.GridFactory()
    >>> factory.insert_vertex([0, 0, 0])
    >>> factory.insert_vertex([1, 0, 0])
    >>> factory.insert_vertex([0, 1, 0])
    >>> factory.insert_element([0, 1, 2])
    >>> grid = factory.finalize()

    """

    def __init__(self):
        """Construct a new GridFactory object."""
        from bempp.core.grid.grid_factory import GridFactory as Factory
        self._impl = Factory()

    def insert_vertex(self, vertex):
        """
        Insert a vertex into a grid.

        The vertex must be a list type object with 3 components.

        """
        self._impl.insert_vertex(vertex)

    def insert_element(self, element, domain_index=0):
        """
        Insert an element into a list.

        The element must be a list type object with three components
        specifying the insertion indices of the three vertices associated
        with the element. The domain_index allows to group elements
        into different sets with different indicies.

        """
        self._impl.insert_element(element, domain_index)

    def finalize(self):
        """Finalize the grid creation and return a grid object."""
        from bempp.api.grid.grid import Grid
        from bempp.api import LOGGER

        grid = Grid(self._impl.finalize())

        LOGGER.info(
            "Created grid with {0} elements, {1} nodes and {2} edges.".format(
                grid.leaf_view.entity_count(0),
                grid.leaf_view.entity_count(2),
                grid.leaf_view.entity_count(1)))
        return grid
