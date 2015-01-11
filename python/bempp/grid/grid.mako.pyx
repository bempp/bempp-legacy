from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp cimport bool as cbool
from libcpp.vector cimport vector
from bempp.utils cimport catch_exception
from bempp.utils cimport unique_ptr
from bempp.utils.armadillo cimport Col, Mat
from bempp.grid.grid_view cimport c_GridView, GridView
from bempp.grid.grid_view cimport _grid_view_from_unique_ptr
import numpy as np
cimport numpy as np
cimport cython


cdef extern from "<array>" namespace "std":
    cdef cppclass c_vertex "std::array<double,3>":
        c_vertex()
        double& operator[](int)

cdef extern from "<array>" namespace "std":
    cdef cppclass c_element "std::array<int,3>":
        c_element()
        int& operator[](int)

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

cdef extern from "bempp/grid/py_sphere.hpp" namespace "Bempp":

    cdef cppclass SphereMesh:
            SphereMesh(int n, double radius, c_vertex offset);
            vector[c_vertex]* nodes() 
            vector[c_element]* elements() 


cdef class Grid:
    """ 
        


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


def grid_from_element_data(vertices, elements, domain_indices=[]):
    """

    Create a grid from a given set of vertices and elements.

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
    A bempp.Grid object.

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
    return grid

def structured_grid(lower_left,upper_right,subdivisions):
    """

    Create a two dimensional grid by defining the lower left and
    upper right point.

    Parameters
    ----------
    lower_left : tuple
        The (x,y) coordinate of the lower left corner of the grid.
    upper_right : tuple
        The (x,y) coordinate of the upper right corner of the grid.
    subdivisions : tuple
        A tuple (N,M) specifiying the number of subdivisions in
        each dimension.

    Returns
    -------
    A bempp.Grid object.

    Examples
    --------
    The following command creates a grid of the unit square [0,1]^2.

    >>> grid = structured_grid((0,0),(1,1),(100,100))

    """


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
    return grid

@cython.boundscheck(False)
@cython.wraparound(False)
def grid_from_sphere(int n, double radius=1.0, object origin = [0,0,0]):
    """

    Create a grid discretizing a sphere.

    Parameters
    ----------
    n : int
        The recursion level for the sphere discretization. The
        grid will have 8*(4**n) elements.
    radius : float
        The radius of the sphere (default 1.0).
    origin : list
        List with 3 elements defining the origin of the sphere
        (default [0,0,0]).

    Returns
    -------
    A bempp.grid object

    Examples
    --------
    To discretize a sphere with radius 2 and
    origin [0,1,0] using 512 elements use

    >>> grid = grid_from_sphere(3,2,[0,1,0])

    """


    cdef c_vertex offset_data;
    cdef SphereMesh* mesh
    cdef vector[c_vertex]* c_nodes
    cdef vector[c_element]* c_elements 
    cdef np.ndarray nodes
    cdef np.ndarray elements
    cdef np.ndarray[double,ndim=2] nodes_buf
    cdef np.ndarray[long,ndim=2] elements_buf

    cdef int number_of_nodes;
    cdef int number_of_elements;

    for i in range(3):
        offset_data[i] = origin[i]

    mesh = new SphereMesh(n,radius,offset_data)
    c_nodes = deref(mesh).nodes()
    c_elements = deref(mesh).elements()

    number_of_nodes = deref(c_nodes).size()
    number_of_elements = deref(c_elements).size()

    nodes = np.empty((3,number_of_nodes),dtype='float64')
    elements = np.empty((3,number_of_elements),dtype='int64')

    nodes_buf = nodes
    elements_buf = elements


    for j in range(3):
        for i in range(number_of_nodes):
            nodes_buf[j,i] = deref(c_nodes)[i][j]
        for i in range(number_of_elements):
            elements_buf[j,i] = deref(c_elements)[i][j]

    del mesh

    return grid_from_element_data(nodes,elements)


        
