#cython: embedsignature=True

from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp cimport bool as cbool
from libcpp.vector cimport vector
from bempp.core.utils cimport catch_exception
from bempp.core.utils cimport unique_ptr
from bempp.core.utils cimport Vector
from bempp.core.utils cimport eigen_matrix_to_np_int
from bempp.core.grid.entity cimport Entity0
from bempp.core.grid.grid_view cimport c_GridView, GridView
from bempp.core.grid.grid_view cimport _grid_view_from_unique_ptr
from bempp.core.grid.id_set cimport IdSet
from bempp.core.cuda cimport CudaGridFloat
from bempp.core.cuda cimport CudaGridDouble
import numpy as _np
cimport numpy as _np
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

    shared_ptr[c_Grid] cart_grid "Bempp::GridFactory::createStructuredGrid"(
            const GridParameters& params,
            const double* lowerLeft,
            const double*  upperRight,
            const int* nElements
    ) except +catch_exception

    shared_ptr[c_Grid] connect_grid \
            "Bempp::GridFactory::createGridFromConnectivityArrays"(
            const GridParameters& params,
            const double* vertices,
            int nvertices,
            const int* elementCorners,
            int nelements,
            vector[int]& domainIndices
    ) except +catch_exception

cdef extern from "bempp/core/grid/py_sphere.hpp" namespace "Bempp":

    cdef cppclass SphereMesh:
            SphereMesh(int n, double radius, c_vertex offset);
            vector[c_vertex]* nodes() 
            vector[c_element]* elements() 


cdef class Grid:
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
    :mod:`bempp.api.shapes`, :func:`bempp.api.import_grid`,
    :class:`bempp.api.GridFactory`
    
    """


    def __cinit__(self):
        self.impl_.reset()

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()

    def __richcmp__(Grid self, Grid other not None, int op):
        if op != 2:
            raise AttributeError("Incorrect operator")
        return self.impl_.get() == other.impl_.get()

    def plot(self):
        """

        Plot a grid with Gmsh.

        """

        from bempp.api.external.viewers import visualize_with_gmsh

        visualize_with_gmsh(self)

    cpdef unsigned int vertex_insertion_index(self, Entity2 vertex):
        return deref(self.impl_).vertexInsertionIndex(
                deref(vertex.impl_))

    cpdef unsigned int element_insertion_index(self, Entity0 element):
        return deref(self.impl_).elementInsertionIndex(
                deref(element.impl_))
    
    cpdef Entity0 element_from_insertion_index(self, int index):

        if self._insertion_index_to_element is None:
            self._insertion_index_to_element = {}
            for element in self.leaf_view.entity_iterator(0):
                index = self.element_insertion_index(element)
                self._insertion_index_to_element[index] = element

        return self._insertion_index_to_element[index]

    cpdef Entity2 vertex_from_insertion_index(self, int index):

        if self._insertion_index_to_vertex is None:
            self._insertion_index_to_vertex = {}
            for vertex in self.leaf_view.entity_iterator(2):
                index = self.vertex_insertion_index(vertex)
                self._insertion_index_to_vertex[index] = vertex 

        return self._insertion_index_to_vertex[index]



    def clone(self):
        """ Create a new grid from a copy of the current grid. """

        vertices = self.leaf_view.vertices
        elements = self.leaf_view.elements
        domain_indices = self.leaf_view.domain_indices
        return grid_from_element_data(vertices, elements, domain_indices)

    def get_mark(self, Entity0 e):
        """ Return the mark status of a given element. """

        return deref(self.impl_).getMark(deref(e.impl_))

    def mark(self, Entity0 e, int refcount = 1):
        """ Mark an element for refinement. """

        deref(self.impl_).mark(refcount, deref(e.impl_))

    def refine(self):
        """ Refine grid. """

        self._grid_view = None
        deref(self.impl_).preAdapt()
        deref(self.impl_).adapt()
        deref(self.impl_).postAdapt()
        deref(self.impl_).sendUpdateSignal()

    def global_refine(self, int refcount = 1):
        """ Refine all elements. """

        self._grid_view = None
        deref(self.impl_).globalRefine(refcount)
        deref(self.impl_).sendUpdateSignal()

    def barycentric_grid(self):
        """Return a barycentrically refined grid."""

        cdef Grid grid = Grid()
        grid.impl_.assign(deref(self.impl_).barycentricGrid())
        return grid

    def barycentric_descendents_map(self):
        """Return the map between elements in the original grid and its barycentric refinement."""

        return eigen_matrix_to_np_int(deref(self.impl_).barycentricSonMap()) 

    def push_to_device(self, device_id, precision='double'):
        """Push a grid to a Cuda device with a given id and precision either 'float' or 'double'."""

        cdef CudaGridFloat cuda_grid_float = CudaGridFloat()
        cdef CudaGridDouble cuda_grid_double = CudaGridDouble()
        if precision == 'double':
            cuda_grid_double.impl_.assign(deref(self.impl_).pushToDeviceDouble(device_id))
            return cuda_grid_double
        elif precision == 'float':
            cuda_grid_float.impl_.assign(deref(self.impl_).pushToDeviceFloat(device_id))
            return cuda_grid_float
        else:
            raise ValueError("precision must be either 'double' or 'float'".)
        return cuda_grid

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
            cdef:
                int n = deref(self.impl_).dimWorld()
                Vector[double] lower
                Vector[double] upper
                int i


            deref(self.impl_).getBoundingBox(lower, upper)
            if upper.rows() != n or lower.rows() != n:
                raise RuntimeError("Error in getBoundingBox")

            cdef _np.ndarray result = _np.ones((2, n), dtype="double", order='C')
            for i in range(n):
                result[0, i] = lower.value(i)
                result[1, i] = upper.value(i)
            return result

    property leaf_view:
        def __get__(self):
            """Return a leaf view onto the Grid"""
            if self._grid_view is not None:
                return self._grid_view

            cdef unique_ptr[c_GridView] view = deref(self.impl_).leafView()
            cdef GridView grid_view = _grid_view_from_unique_ptr(view)
            grid_view._grid = self
            self._grid_view = grid_view
            return grid_view

    property id_set:
        def __get__(self):
            """Return an Id set for the grid"""

            cdef IdSet id_set = IdSet()
            id_set.impl_ = &(deref(self.impl_).globalIdSet())
            id_set._grid = self
            return id_set


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


    from numpy import require
    cdef:
        GridParameters parameters
        double[::1, :] vert_ptr \
                = require(vertices, "double", 'F')
        int[::1, :] corners_ptr = require(elements, "intc", 'F')
        vector[int] indices
        Grid grid = Grid.__new__(Grid)
    for index in domain_indices:
        indices.push_back(int(index))
    parameters.topology = TRIANGULAR
    try:
        grid.impl_ = connect_grid(parameters, &vert_ptr[0,0],
                vert_ptr.shape[1], &corners_ptr[0,0], 
                corners_ptr.shape[1], indices)
    except:
        raise
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
    grid : bempp.Grid
        A structured grid.

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
    if len(ll_ptr) != 2 or len(ur_ptr) != 2 or len(n_ptr)!= 2:
        raise ValueError("Inputs have wrong dimension.")

    cdef Grid grid = Grid.__new__(Grid)
    parameters.topology = TRIANGULAR
    try:
        grid.impl_ = cart_grid(parameters, &ll_ptr[0], &ur_ptr[0],
                 &n_ptr[0])
    except:
        raise
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
    grid : bempp.Grid
        The discretization of a sphere.

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
    cdef _np.ndarray nodes
    cdef _np.ndarray elements
    cdef _np.ndarray[double,ndim=2] nodes_buf
    cdef _np.ndarray[int,ndim=2] elements_buf

    cdef int number_of_nodes;
    cdef int number_of_elements;

    for i in range(3):
        offset_data[i] = origin[i]

    mesh = new SphereMesh(n,radius,offset_data)
    c_nodes = deref(mesh).nodes()
    c_elements = deref(mesh).elements()

    number_of_nodes = deref(c_nodes).size()
    number_of_elements = deref(c_elements).size()

    nodes = _np.empty((3,number_of_nodes),dtype='float64')
    elements = _np.empty((3,number_of_elements),dtype='intc')

    nodes_buf = nodes
    elements_buf = elements


    for j in range(3):
        for i in range(number_of_nodes):
            nodes_buf[j,i] = deref(c_nodes)[i][j]
        for i in range(number_of_elements):
            elements_buf[j,i] = deref(c_elements)[i][j]

    del mesh

    return grid_from_element_data(nodes,elements)





        
