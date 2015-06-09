from bempp.utils cimport shared_ptr,unique_ptr
from bempp.utils cimport Vector
from bempp.grid.grid_view cimport c_GridView, GridView
from bempp.grid.entity cimport Entity0, Entity2
from bempp.grid.entity cimport c_Entity
from bempp.grid.codim_template cimport codim_zero,codim_one,codim_two

cdef extern from "bempp/grid/grid.hpp" namespace "Bempp" nogil:
    cdef cppclass c_Grid "Bempp::Grid":
        int dim() const
        int dimWorld() const
        int maxLevel() const
        int topology() const
        unique_ptr[c_GridView] leafView() const
        void getBoundingBox(const Vector[double]&, const Vector[double]&) const
        unsigned int vertexInsertionIndex(const c_Entity[codim_two]&) const
        unsigned int elementInsertionIndex(const c_Entity[codim_zero]&) const

    cdef enum Topology "Bempp::GridParameters::Topology":
        LINEAR "Bempp::GridParameters::LINEAR"
        TRIANGULAR "Bempp::GridParameters::TRIANGULAR"
        QUADRILATERAL "Bempp::GridParameters::QUADRILATERAL"
        HYBRID_2D "Bempp::GridParameters::HYBRID_2D"
        TETRAHEDRAL "Bempp::GridParameters::TETRAHEDRAL"

    cdef cppclass GridParameters:
        Topology topology

cdef class Grid:
    ## Holds pointer to C++ implementation
    cdef shared_ptr[const c_Grid] impl_
    cdef GridView _grid_view
    cpdef unsigned int vertex_insertion_index(self,Entity2 vertex)
    cpdef unsigned int element_insertion_index(self,Entity0 element)



cdef Grid mark_and_refine(Grid grid, object index_list)
