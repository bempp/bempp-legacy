from bempp.core.utils cimport unique_ptr
from bempp.core.utils cimport catch_exception
from cython.operator cimport dereference as deref
from cython cimport address
from libcpp cimport bool as cbool

from bempp.core.grid.entity cimport c_Entity
from bempp.core.grid.entity_pointer cimport c_EntityPointer

from bempp.core.grid.entity cimport Entity0
from bempp.core.grid.entity cimport Entity1
from bempp.core.grid.entity cimport Entity2
from bempp.core.grid.entity_pointer cimport EntityPointer0
from bempp.core.grid.entity_pointer cimport EntityPointer1
from bempp.core.grid.entity_pointer cimport EntityPointer2

cdef class EntityIterator0:


    def __dealloc__(self):
        self.impl_.reset()

    cdef cbool finished(self):
        return deref(self.impl_).finished()

    cdef void c_next(self) except + :
        deref(self.impl_).next()

    cdef EntityPointer0 _frozen(self):
        cdef EntityPointer0 ep = EntityPointer0()
        cdef unique_ptr[c_EntityPointer[codim_zero]] c_ep = deref(self.impl_).frozen()
        ep.impl_.swap(c_ep)
        return ep

    def __next__(self):
        cdef EntityPointer0 ep
        if not self.finished():
            ep = self._frozen()
            self.c_next()
            return ep.entity
        else:
            raise StopIteration()

    def __iter__(self):
        return self

cdef class EntityIterator1:

    def __dealloc__(self):
        self.impl_.reset()

    cdef cbool finished(self):
        return deref(self.impl_).finished()

    cdef void c_next(self) except + :
        deref(self.impl_).next()

    cdef EntityPointer1 _frozen(self):
        cdef EntityPointer1 ep = EntityPointer1()
        cdef unique_ptr[c_EntityPointer[codim_one]] c_ep = deref(self.impl_).frozen()
        ep.impl_.swap(c_ep)
        return ep

    def __next__(self):
        cdef EntityPointer1 ep
        if not self.finished():
            ep = self._frozen()
            self.c_next()
            return ep.entity
        else:
            raise StopIteration()

    def __iter__(self):
        return self


cdef class EntityIterator2:

    def __dealloc__(self):
        self.impl_.reset()

    cdef cbool finished(self):
        return deref(self.impl_).finished()

    cdef void c_next(self) except + :
        deref(self.impl_).next()

    cdef EntityPointer2 _frozen(self):
        cdef EntityPointer2 ep = EntityPointer2()
        cdef unique_ptr[c_EntityPointer[codim_two]] c_ep = deref(self.impl_).frozen()
        ep.impl_.swap(c_ep)
        return ep

    def __next__(self):
        cdef EntityPointer2 ep
        if not self.finished():
            ep = self._frozen()
            self.c_next()
            return ep.entity
        else:
            raise StopIteration()

    def __iter__(self):
        return self



       




        


