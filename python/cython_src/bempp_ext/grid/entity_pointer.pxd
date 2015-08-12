from bempp_ext.utils cimport unique_ptr
from bempp_ext.grid.codim_template cimport codim_zero,codim_one,codim_two

from bempp_ext.grid.entity cimport c_Entity

cdef extern from "bempp/grid/entity_pointer.hpp" namespace "Bempp":
    cdef cppclass c_EntityPointer "Bempp::EntityPointer"[codim]:
         const int codimension
         const c_Entity[codim]& entity() const


from bempp_ext.grid.entity cimport Entity0
from bempp_ext.grid.entity cimport Entity1
from bempp_ext.grid.entity cimport Entity2

cdef class EntityPointer0:
    cdef unique_ptr[c_EntityPointer[codim_zero]] impl_
    cpdef Entity0 _entity(self)

cdef class EntityPointer1:
    cdef unique_ptr[c_EntityPointer[codim_one]] impl_
    cpdef Entity1 _entity(self)

cdef class EntityPointer2:
    cdef unique_ptr[c_EntityPointer[codim_two]] impl_
    cpdef Entity2 _entity(self)
