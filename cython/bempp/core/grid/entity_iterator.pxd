cimport numpy
from bempp.core.utils cimport unique_ptr, catch_exception
from libcpp cimport bool as cbool
from bempp.core.grid.codim_template cimport codim_zero,codim_one,codim_two
from bempp.core.grid.entity_pointer cimport c_EntityPointer

cdef extern from "bempp/grid/entity_iterator.hpp" namespace "Bempp":
    cdef cppclass c_EntityIterator "Bempp::EntityIterator"[codim]:
        void next()
        cbool finished() const
        unique_ptr[c_EntityPointer[codim]] frozen() const

from bempp.core.grid.entity_pointer cimport EntityPointer0
from bempp.core.grid.entity_pointer cimport EntityPointer1
from bempp.core.grid.entity_pointer cimport EntityPointer2

cdef class EntityIterator0:
    cdef unique_ptr[c_EntityIterator[codim_zero]] impl_
    cdef cbool finished(self)
    cdef void c_next(self) except +
    cdef EntityPointer0 _frozen(self)

cdef class EntityIterator1:
    cdef unique_ptr[c_EntityIterator[codim_one]] impl_
    cdef cbool finished(self)
    cdef void c_next(self) except +
    cdef EntityPointer1 _frozen(self)

cdef class EntityIterator2:
    cdef unique_ptr[c_EntityIterator[codim_two]] impl_
    cdef cbool finished(self)
    cdef void c_next(self) except +
    cdef EntityPointer2 _frozen(self)
