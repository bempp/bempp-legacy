"""
Definition of the Grid class and factory functions.

The grid stores all relevant information to query grids and
iterate over entities. BEM++ Grids use the Dune-Grid library for
their implementation. For a more detailed description of the
underlying concepts see therefore also the tutorials at
https://dune-project.original

"""
from bempp.api.grid.entity import Entity as _Entity
from bempp.api.grid.grid_view import GridView as _GridView
from bempp.api.grid.id_set import IdSet as _IdSet

#pylint: disable=protected-access

class Grid(object):
    """The basic grid class in BEM++.

    The Grid class contains interfaces to query and iterate grid
    information in BEM++. The basic access point to iterate
    through grid data structures is the method `leaf_view` that
    returns a GridView object. The BEM++ Grid class is implemented
    using using the Dune-grid library and follows closely its
    design principles.

    The Grid class can not be directly instantiated. Rather,
    a range of factory functions are available that return grids.

    See Also
    --------
    :mod:`bempp.api.shapes`, :func:`bempp.api.grid.import_grid`,
    :class:`bempp.api.grid.GridFactory`

    """

    def __init__(self, impl):
        """Constructor. Not to be called diretly."""
        self._impl = impl
        self._barycentric_grid = None
        self._view = None
        self._marked_set = {}

    def plot(self):
        """Plot the grid with Gmsh."""
        from bempp.api.external.viewers import visualize
        visualize(self)

    def vertex_insertion_index(self, vertex):
        """Return the vertex insertion index of the given vertex."""
        return self._impl.vertex_insertion_index(vertex._impl)

    def element_insertion_index(self, element):
        """Return the element insertion index of the given element."""
        return self._impl.element_insertion_index(element._impl)

    def element_from_insertion_index(self, index):
        """Return the element associated with a given insertion index."""
        return _Entity(0, self._impl.element_from_insertion_index(index))

    def vertex_from_insertion_index(self, index):
        """Return the vertex associated with a given insertion index."""
        return _Entity(2, self._impl.vertex_from_insertion_index(index))

    def clone(self):
        """Create a new grid from a copy of the current grid."""
        return Grid(self._impl.clone())

    def get_mark(self, element):
        """Return the mark status of a given element."""
        return self.leaf_view.get_mark(element._impl)

    def mark(self, element):
        """Mark an element for refinement."""
        self.leaf_view.mark(element)

    def refine(self):
        """Refine grid."""
        return self.leaf_view.refine()

    def clear_marks(self):
        """Clear all marks."""
        self.leaf_view.clear_marks()

    def barycentric_grid(self):
        """Return a barycentrically refine grid."""
        if self._barycentric_grid is None:
            self._barycentric_grid = Grid(self._impl.barycentric_grid())
        return self._barycentric_grid

    def barycentric_descendents_map(self):
        """
        Index map between a grid and its barycentric refinement.

        Return a matrix that provides a map between elements in the original
        grid and the barycentric refinement.

        This function returns a (nelements x 6) matrix where the row i contains
        the indices of all elements in the barycentric grid that are
        descendents of the element with index i in the original grid.

        """
        return self._impl.barycentric_descendents_map()

    @property
    def dim(self):
        """Return the dimension of the grid."""
        return self._impl.dim

    @property
    def dim_world(self):
        """Return the dimension of the space containing the grid."""
        return self._impl.dim_world

    @property
    def max_level(self):
        """Return the maximum level defind in this grid."""
        return self._impl.max_level

    @property
    def topology(self):
        """Return the Grid topology."""

    @property
    def bounding_box(self):
        """Return a (2 x 3) array with the bounding box of the grid.

        Examples
        --------
        If the grid is a unit cube the bounding box is the array

        [[0, 0, 0]
         [1, 1, 1]].

        The columns are the dimensions x, y and z. The first row
        is the lower bound and the second row is the upper bound.

        """
        return self._impl.bounding_box

    @property
    def leaf_view(self):
        """Return a view onto the grid."""
        if self._view is None:
            self._view = _GridView(self._impl.leaf_view)
        return self._view

    @property
    def id_set(self):
        """Return an Id set for the grid."""
        return _IdSet(self._impl.id_set)


def grid_from_element_data(vertices, elements, domain_indices=None):
    """Create a grid from a given set of vertices and elements.

    This function takes a list of vertices and a list of elements
    and returns a grid object.

    Parameters
    ----------
    vertices : np.ndarray[float]
        A (3xN) array of vertices.
    elements : np.ndarray[int]
        A (3xN) array of elements.

    Returns
    -------
    grid : bempp.Grid
        The grid representing the specified element data.

    Examples
    --------
    The following code creates a grid with two elements.

    >>> import numpy as np
    >>> vertices = np.array([[0,1,1,0],
                             [0,0,1,1],
                             [0,0,0,0]])
    >>> elements = np.array([[0,1],
                             [1,2],
                             [3,3]])
    >>> grid = grid_from_element_data(vertices,elements)

    """
    from bempp.api import log
    #pylint: disable=no-name-in-module
    from bempp.core.grid.grid import grid_from_element_data as grid_fun

    if domain_indices is None:
        domain_indices = []

    grid = Grid(grid_fun(vertices, elements, domain_indices))
    log(
        "Created grid with %i elements, %i nodes and %i edges." %(
        grid.leaf_view.entity_count(0),
        grid.leaf_view.entity_count(2),
        grid.leaf_view.entity_count(1)))
    return grid


def structured_grid(lower_left, upper_right,
                    subdivisions, axis="xy", offset=0):
    """
    Create a two dimensional rectangular grid.

    Parameters
    ----------
    lower_left : tuple
        The (x,y) coordinate of the lower left corner of the grid.
    upper_right : tuple
        The (x,y) coordinate of the upper right corner of the grid.
    subdivisions : tuple
        A tuple (N,M) specifiying the number of subdivisions in
        each dimension.
    axis : string
        Possible choices are "xy", "xz", "yz". Denotes the
        axes along which the structured grid is generated.
        Default is "xy".
    offset : double
        Defines an offset value that shifts the structured grid
        in the remaining coordinate direction.

    Returns
    -------
    grid : bempp.Grid
        A structured grid.

    Examples
    --------
    The following command creates a grid of the unit square [0,1]^2.

    >>> grid = structured_grid((0,0),(1,1),(100,100))

    """
    import numpy as np
    #pylint: disable=no-name-in-module
    from bempp.core.grid.grid import structured_grid as grid_fun

    # Get a grid along the xy axis
    grid = Grid(grid_fun(lower_left, upper_right, subdivisions))

    vertices = grid.leaf_view.vertices
    elements = grid.leaf_view.elements

    #pylint: disable=no-member

    if axis == "xy":
        # Nothing to be done
        offset_vector = np.array([[0, 0, offset]]).T
        return grid_from_element_data(vertices + offset_vector, elements)
    elif axis == "xz":
        offset_vector = np.array([[0, offset, 0]]).T
        return grid_from_element_data(
            vertices[[0, 2, 1], :] + offset_vector, elements)
    elif axis == "yz":
        offset_vector = np.array([[offset, 0, 0]]).T
        return grid_from_element_data(
            vertices[[2, 0, 1], :] + offset_vector, elements)
    else:
        raise ValueError("axis parameter must be 'xy', 'xz', or 'yz'.")
