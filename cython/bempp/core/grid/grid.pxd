from bempp.core.utils cimport shared_ptr,unique_ptr
from bempp.core.utils cimport Vector
from bempp.core.utils cimport Matrix
from bempp.core.grid.grid_view cimport c_GridView, GridView
from bempp.core.grid.entity cimport Entity0, Entity2
from bempp.core.grid.entity cimport c_Entity
from bempp.core.grid.codim_template cimport codim_zero,codim_one,codim_two
from bempp.core.grid.id_set cimport c_IdSet
from libcpp cimport bool as cbool

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
        cbool mark(int refCount, const c_Entity[codim_zero]&)
        cbool adapt()
        cbool preAdapt()
        void postAdapt()
        void sendUpdateSignal() const
        void globalRefine(int refCount)
        int getMark(const c_Entity[codim_zero]&)
        c_IdSet& globalIdSet() const
        shared_ptr[c_Grid] barycentricGrid() const
        Matrix[int] barycentricSonMap() const


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
    cdef shared_ptr[c_Grid] impl_
    cdef GridView _grid_view
    cdef object _insertion_index_to_element
    cdef object _insertion_index_to_vertex
    cpdef unsigned int vertex_insertion_index(self,Entity2 vertex)
    cpdef unsigned int element_insertion_index(self,Entity0 element)
    cpdef Entity0 element_from_insertion_index(self, int index)
    cpdef Entity2 vertex_from_insertion_index(self, int index)


