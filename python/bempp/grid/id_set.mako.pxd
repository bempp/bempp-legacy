from bempp.grid.codim_template cimport codim_zero,codim_one,codim_two
from bempp.grid.entity cimport c_Entity
from bempp.grid.entity cimport Entity0, Entity1, Entity2

cdef extern from "bempp/grid/id_set.hpp" namespace "Bempp":
    cdef cppclass c_IdSet "Bempp::IdSet":
        size_t entityId(const c_Entity[codim_zero]&) const
        size_t entityId(const c_Entity[codim_one]&) const
        size_t entityId(const c_Entity[codim_two]&) const
        size_t subEntityId(const c_Entity[codim_zero]&, size_t i,
                int codimSub) const

from bempp.grid.grid cimport Grid

cdef class IdSet:
    cdef Grid _grid
    cdef const c_IdSet* impl_

    cdef size_t entity_id_0(self, Entity0 entity)
    cdef size_t entity_id_1(self, Entity1 entity)
    cdef size_t entity_id_2(self, Entity2 entity)

    cpdef size_t sub_entity_id(self, Entity0 element, size_t i, int codim_sub)




