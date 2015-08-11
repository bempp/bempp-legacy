<%
codims = [('0','codim_zero'),('1','codim_one'),('2','codim_two')]
%>

from bempp.utils cimport unique_ptr
from bempp.grid.codim_template cimport codim_zero,codim_one,codim_two

from bempp.grid.entity cimport c_Entity

cdef extern from "bempp/grid/entity_pointer.hpp" namespace "Bempp":
    cdef cppclass c_EntityPointer "Bempp::EntityPointer"[codim]:
         const int codimension
         const c_Entity[codim]& entity() const


% for (codim,codim_template) in codims:
from bempp.grid.entity cimport Entity${codim}

cdef class EntityPointer${codim}:
    cdef unique_ptr[c_EntityPointer${'['+codim_template+']]'} impl_
    cpdef Entity${codim} _entity(self)

% endfor
