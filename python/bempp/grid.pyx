from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp cimport bool as cbool
from bempp.utils cimport catch_exception

cdef extern from "bempp/grid/grid.hpp" namespace "Bempp":
    cdef enum Topology "Bempp::GridParameters::Topology":
        LINEAR "Bempp::GridParameters::LINEAR"
        TRIANGULAR "Bempp::GridParameters::TRIANGULAR"
        QUADRILATERAL "Bempp::GridParameters::QUADRILATERAL"
        HYBRID_2D "Bempp::GridParameters::HYBRID_2D"
        TETRAHEDRAL "Bempp::GridParameters::TRIANGULAR"

    cdef cppclass GridParameters:
        Topology topology


cdef extern from "bempp/grid/grid_factory.hpp" namespace "Bempp":
    shared_ptr[c_Grid] grid_from_file "Bempp::GridFactory::importGmshGrid"(
            const GridParameters& params,
            const string& fileName,
            cbool verbose,
            cbool insertBoundarySegments
    ) except +catch_exception


cdef class Grid:
    def __cinit__(self, **kwargs):
        cdef:
            GridParameters grid_parameters
            string filename
        from_file = kwargs.get("topology", None) is not None \
            and kwargs.get("filename", None) is not None
        if from_file:
            topology = {
                'linear': LINEAR, 'l': LINEAR,
                'triangular': TRIANGULAR, 'triangle': TRIANGULAR,
                'quadrilateral': QUADRILATERAL, 'q': QUADRILATERAL,
                'hybrid2d': HYBRID_2D, 'hybrid_2d': HYBRID_2D, 'h': HYBRID_2D,
                'tetrahedral': TETRAHEDRAL, 'tetra': TETRAHEDRAL
            }[kwargs.pop('topology').lower()]
            grid_parameters.topology = topology
            filename = str(kwargs.pop('filename'))
            self.impl_ = grid_from_file(
                    grid_parameters, filename,
                    kwargs.pop('verbose', False) == True,
                    kwargs.pop('insert_boundary_segments', False) == True
            )
        else: raise RuntimeError("Incorrect set of arguments to Grid")

    property dim:
        def __get__(self):
            return deref(self.impl_).dim()
    property dim_world:
        def __get__(self):
            return deref(self.impl_).dimWorld()
    property max_level:
        def __get__(self):
            return deref(self.impl_).maxLevel()
    property topology:
        def __get__(self):
            cdef int value = deref(self.impl_).topology()
            if value == LINEAR: return 'linear'
            elif value == TRIANGULAR: return 'triangular'
            elif value == QUADRILATERAL: return 'quadrilateral'
            elif value == HYBRID_2D: return 'hybrid2d'
            elif value == TETRAHEDRAL: return 'tetrahedral'
            raise RuntimeError("C++ to Python bug: Unknown topology")
