from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp cimport bool as cbool
from bempp.utils cimport catch_exception
from bempp.utils.armadillo cimport Col

cdef extern from "bempp/grid/grid_factory.hpp" namespace "Bempp":
    shared_ptr[c_Grid] grid_from_file "Bempp::GridFactory::importGmshGrid"(
            GridParameters& params,
            string& fileName,
            cbool verbose,
            cbool insertBoundarySegments
    ) except +catch_exception

    shared_ptr[c_Grid] cart_grid "Bempp::GridFactory::createStructuredGrid"(
            GridParameters& params,
            Col[double]& lowerLeft,
            Col[double]& upperRight,
            Col[unsigned int]& nElements
    ) except +catch_exception


cdef class Grid:
    """ Grid information

        Initialization can be done via several sets of parameters:

        1. From a file:
            - topology
            - filename
        1. Structure grid:
            - topology
            - lower_left
            - upper_right
            - nelements

        Parameters:
        -----------

        topology : string
            The topology of the grid. One of "linear", "triangular"
            "quadrilateral", "hybrid_2d", or "tetrahedral".

        filename : str
            Path to the file. The path is expanded for environment
            variables and such.

        upper_right : tuple
            Coordinates of the upper right coordinate of a structured grid.

        lower_left : tuple
            Coordinates of the lower left coordinate of a structured grid.

        nelements : tuple
            Number of elements for each axis of a structured grid.

    """
    cdef void __create_from_file(self, dict kwargs,
            GridParameters &parameters) except *:
        from os.path import abspath, expandvars, expanduser
        cdef string cfilename
        cfilename = abspath(expandvars(expanduser(kwargs['filename'])))
        self.impl_ = grid_from_file(
                parameters, cfilename,
                kwargs.pop('verbose', False) == True,
                kwargs.pop('insert_boundary_segments', False) == True
        )

    cdef void __create_cartesian_grid(self, dict kwargs,
            GridParameters& parameters) except *:
        from numpy import require
        lower_left = require(kwargs['lower_left'], "double", 'C')
        upper_right = require(kwargs['upper_right'], "double", 'C')
        subdivisions = require(kwargs['subdivisions'], "intc", 'C')
        nelements = len(lower_left)
        for o in [lower_left, upper_right, subdivisions]:
            if not len(o) == nelements:
                raise ValueError("Inputs have differing lengths")
            if len(o.shape) != 1:
                raise ValueError("Inputs are not vectors")
        cdef:
            double[::1] ll_ptr = lower_left
            double[::1] ur_ptr = upper_right
            int[::1] n_ptr = subdivisions
            Col[double]* c_lower_left = new Col[double](&ll_ptr[0],
                nelements, False, True)
            Col[double]* c_upper_right = new Col[double](&ur_ptr[0],
                nelements, False, True)
            Col[unsigned int]* c_subdivisions = \
                    new Col[unsigned int](<unsigned int*>&n_ptr[0],
                            nelements, False, True)
        try:
            self.impl_ = cart_grid(parameters, deref(c_lower_left),
                    deref(c_upper_right), deref(c_subdivisions))
        except:
            del c_lower_left
            del c_upper_right
            del c_subdivisions
            raise

    def __cinit__(self, topology='triangular', **kwargs):
        cdef GridParameters parameters

        try:
            parameters.topology = {
                'linear': LINEAR, 'l': LINEAR,
                'triangular': TRIANGULAR, 'triangle': TRIANGULAR,
                'quadrilateral': QUADRILATERAL, 'q': QUADRILATERAL,
                'hybrid2d': HYBRID_2D, 'hybrid_2d': HYBRID_2D, 'h': HYBRID_2D,
                'tetrahedral': TETRAHEDRAL, 'tetra': TETRAHEDRAL
            }[topology]
        except KeyError:
            raise ValueError("Incorrect topology %s" % topology)

        # Check which set of input has been given
        check = lambda x: all([kwargs.get(u, None) for u in x])
        use_file = check(['filename'])
        use_structured = check(['lower_left', 'upper_right', 'subdivisions'])

        # Makes sure at least one set of argument is provided
        nargs = sum([1 if x else 0 for x in [use_file, use_structured]])
        if nargs == 0:
            raise TypeError("Incorrect set of input arguments")
        elif nargs > 1:
            msg = "Ambiguous input: matches %i set of input arguments"
            raise TypeError(msg % nargs)

        # At this point, we can choose how to initialize the grid
        if use_file:
            self.__create_from_file(kwargs, parameters)
        elif use_structured:
            self.__create_cartesian_grid(kwargs, parameters)



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
            if value == LINEAR: return 'linear'
            elif value == TRIANGULAR: return 'triangular'
            elif value == QUADRILATERAL: return 'quadrilateral'
            elif value == HYBRID_2D: return 'hybrid2d'
            elif value == TETRAHEDRAL: return 'tetrahedral'
            raise RuntimeError("C++ to Python bug: Unknown topology")
