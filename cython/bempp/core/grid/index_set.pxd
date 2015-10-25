from bempp.core.grid.codim_template cimport codim_zero,codim_one,codim_two
from bempp.core.grid.entity cimport c_Entity
from bempp.core.grid.entity cimport Entity0, Entity1, Entity2
from bempp.core.grid.grid_view cimport GridView

cdef extern from "bempp/grid/index_set.hpp" namespace "Bempp":
    cdef cppclass c_IndexSet "Bempp::IndexSet":
        size_t entityIndex(const c_Entity[codim_zero]&) const
        size_t entityIndex(const c_Entity[codim_one]&) const
        size_t entityIndex(const c_Entity[codim_two]&) const
        size_t subEntityIndex(const c_Entity[codim_zero]&, size_t i,
                int codimSub) const


cdef class IndexSet:
    cdef GridView _grid_view
    cdef const c_IndexSet* impl_

    cdef size_t entity_index_0(self, Entity0 entity)
    cdef size_t entity_index_1(self, Entity1 entity)
    cdef size_t entity_index_2(self, Entity2 entity)

    cpdef size_t sub_entity_index(self, Entity0 element, size_t i, int codim_sub)




