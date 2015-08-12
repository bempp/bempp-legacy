from bempp_ext.utils cimport unique_ptr
from bempp_ext.grid.geometry cimport c_Geometry
from bempp_ext.grid.codim_template cimport codim_zero,codim_one,codim_two
from libcpp cimport bool as cbool

cdef extern from "bempp/grid/entity.hpp" namespace "Bempp":
    cdef cppclass c_Entity "Bempp::Entity"[codim]:
        size_t level() const
        const c_Geometry& geometry() const
        int domain() # Note: only exists for codim=0
        cbool isLeaf() const # Note: only exists for codim=0
        cbool hasFather() const # Note: only exists for codim=0
        

from bempp_ext.grid.entity_pointer cimport EntityPointer0
from bempp_ext.grid.entity_pointer cimport EntityPointer1
from bempp_ext.grid.entity_pointer cimport EntityPointer2

from bempp_ext.grid.geometry cimport Geometry0
from bempp_ext.grid.geometry cimport Geometry1
from bempp_ext.grid.geometry cimport Geometry2

cdef class Entity0:
   cdef EntityPointer0 _entity_pointer
   cdef const c_Entity[codim_zero] * impl_
   cpdef size_t level(self)
   cpdef Geometry0 _geometry(self)

cdef class Entity1:
   cdef EntityPointer1 _entity_pointer
   cdef const c_Entity[codim_one] * impl_
   cpdef size_t level(self)
   cpdef Geometry1 _geometry(self)

cdef class Entity2:
   cdef EntityPointer2 _entity_pointer
   cdef const c_Entity[codim_two] * impl_
   cpdef size_t level(self)
   cpdef Geometry2 _geometry(self)
