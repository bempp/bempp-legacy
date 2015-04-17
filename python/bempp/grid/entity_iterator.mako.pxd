<%
codims = [('0','codim_zero'),('1','codim_one'),('2','codim_two')]
%>

from bempp.utils cimport unique_ptr, catch_exception
from libcpp cimport bool as cbool
from bempp.grid.codim_template cimport codim_zero,codim_one,codim_two
from bempp.grid.entity_pointer cimport c_EntityPointer

cdef extern from "bempp/grid/entity_iterator.hpp" namespace "Bempp":
    cdef cppclass c_EntityIterator "Bempp::EntityIterator"[codim]:
        void next()
        cbool finished() const
        unique_ptr[c_EntityPointer[codim]] frozen() const

% for (codim,codim_template) in codims:

from bempp.grid.entity_pointer cimport EntityPointer${codim}

cdef class EntityIterator${codim}:
    cdef unique_ptr[c_EntityIterator[${codim_template}]] impl_
    cdef cbool finished(self)
    cdef void c_next(self) except +
    cdef EntityPointer${codim} _frozen(self)

% endfor


