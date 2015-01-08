from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp cimport bool as cbool
from libcpp.vector cimport vector
from bempp.utils cimport catch_exception
from bempp.utils cimport unique_ptr
from bempp.utils.armadillo cimport Col, Mat
from bempp.grid.grid_view cimport c_GridView, GridView
from bempp.grid.grid_view cimport _grid_view_from_unique_ptr

cdef extern from "bempp/grid/grid_factory.hpp" namespace "Bempp":

    shared_ptr[const c_Grid] cart_grid "Bempp::GridFactory::createStructuredGrid"(
            const GridParameters& params,
            Col[double]& lowerLeft,
            Col[double]& upperRight,
            Col[unsigned int]& nElements
    ) except +catch_exception

    shared_ptr[const c_Grid] connect_grid \
            "Bempp::GridFactory::createGridFromConnectivityArrays"(
            const GridParameters& params,
            Mat[double]& vertices,
            Mat[int]& elementCorners,
            vector[int]& domainIndices
    ) except +catch_exception


cdef class Grid:
    """ 
        
Direct initialization can be done via two sets of parameters:

1. Structure grid:
    - lower_left
    - upper_right
    - nelements
2. From numpy arrays:
    - vertices
    - elements
    - domain_indices

To import a grid from a Gmsh file see bempp.file_interfaces.gmsh.


Parameters:
-----------
lower_left : tuple
    Coordinates of the lower left coordinate of a structured grid.

upper_right : tuple
    Coordinates of the upper right coordinate of a structured grid.

nelements : tuple
    Number of elements for each axis of a structured grid.

vertices : numpy.ndarray
    (3xN) Array of vertices.

elements : numpy.ndarray
    (3xN) Array of vertex indices defining the elements.

domain_indicies : list[int]
    Integer list containing the domain indicies of the elements (optional).
    

Examples
--------

To discretise the screen [0,1]x[0,1]x0 use

>>> grid = Grid(lower_left=(0,0),upper_right=(1,1),nelements=(100,100))

The following creates a grid consisting of two elements

>>> grid = Grid(



    """


    def __cinit__(self):
        self.impl_.reset(<const c_Grid*>NULL)

    def __init__(self):
        pass

    def __richcmp__(Grid self, Grid other not None, int op):
        if op != 2:
            raise AttributeError("Incorrect operator")
        return self.impl_.get() == other.impl_.get()


    property dim:
        """" Dimension of the grid. """
        def __get__(self):
            return deref(self.impl_).dim()

    property dim_world:
        """ Dimension of the space containing the grid. """
        def __get__(self):
            return deref(self.impl_).dimWorld()

    property max_level:
        """ Maximum level defined in this grid.

            0 is the coarsest level.
        """
        def __get__(self):
            return deref(self.impl_).maxLevel()

    property topology:
        """ Grid topology """
        def __get__(self):
            cdef int value = deref(self.impl_).topology()
            if value == LINEAR:          return 'linear'
            elif value == TRIANGULAR:    return 'triangular'
            elif value == QUADRILATERAL: return 'quadrilateral'
            elif value == HYBRID_2D:     return 'hybrid2d'
            elif value == TETRAHEDRAL:   return 'tetrahedral'
            raise RuntimeError("C++ to Python bug: Unknown topology")

    property bounding_box:
        """ Bounding box surrounding the grid """
        def __get__(self):
            from numpy import ones
            cdef:
                int n = deref(self.impl_).dimWorld()
                Col[double] lower
                Col[double] upper


            deref(self.impl_).getBoundingBox(lower, upper)
            if upper.n_rows != n or lower.n_rows != n:
                raise RuntimeError("Error in getBoundingBox")

            result = ones((2, n), dtype="double", order='C')
            for i in range(n):
                result[0, i] = lower.at(i)
                result[1, i] = upper.at(i)
            return result

    property leaf_view:
        def __get__(self):
            """Return a leaf view onto the Grid"""
            cdef unique_ptr[c_GridView] view = deref(self.impl_).leafView()
            return _grid_view_from_unique_ptr(view)


def create_grid_from_element_data(vertices, elements, domain_indices=[]):
    from numpy import require
    cdef:
        GridParameters parameters
        double[::1, :] vert_ptr \
                = require(vertices, "double", 'F')
        int[::1, :] corners_ptr = require(elements, "intc", 'F')
        vector[int] indices
        Mat[double]* c_vertices = new Mat[double](&vert_ptr[0, 0],
            vert_ptr.shape[0], vert_ptr.shape[1], False, True)
        Mat[int]* c_corners = new Mat[int](&corners_ptr[0, 0],
            corners_ptr.shape[0], corners_ptr.shape[1], False, True)
        Grid grid = Grid.__new__(Grid)
    for index in domain_indices:
        indices.push_back(int(index))
    parameters.topology = TRIANGULAR
    try:
        grid.impl_ = connect_grid(parameters, deref(c_vertices),
                deref(c_corners), indices)
    except:
        del c_vertices
        del c_corners
        raise
    del c_vertices
    del c_corners

def create_structured_grid(lower_left,upper_right,subdivisions):
    from numpy import require
    cdef:
        GridParameters parameters
        double[::1] ll_ptr = require(lower_left, "double", 'C')
        double[::1] ur_ptr = require(upper_right, "double", 'C')
        int[::1] n_ptr = require(subdivisions, "intc", 'C')
        int nelements = len(ll_ptr)
    if len(ll_ptr) != len(ur_ptr) or len(ll_ptr) != len(n_ptr):
        raise ValueError("Inputs have differing lengths")
    cdef:
        Col[double]* c_lower_left = new Col[double](&ll_ptr[0],
            nelements, False, True)
        Col[double]* c_upper_right = new Col[double](&ur_ptr[0],
            nelements, False, True)
        Col[unsigned int]* c_subdivisions = \
                new Col[unsigned int](<unsigned int*>&n_ptr[0],
                        nelements, False, True)
        Grid grid = Grid.__new__(Grid)
    parameters.topology = TRIANGULAR
    try:
        grid.impl_ = cart_grid(parameters, deref(c_lower_left),
                deref(c_upper_right), deref(c_subdivisions))
    except:
        del c_lower_left
        del c_upper_right
        del c_subdivisions
        raise
    del c_lower_left
    del c_upper_right
    del c_subdivisions

