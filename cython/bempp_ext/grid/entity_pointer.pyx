from cython cimport address
from cython.operator cimport dereference as deref

from bempp_ext.grid.entity cimport Entity0
from bempp_ext.grid.entity cimport Entity1
from bempp_ext.grid.entity cimport Entity2

cdef class EntityPointer0:

    cpdef Entity0 _entity(self):
        cpdef Entity0 e = Entity0()
        e.impl_ = address(deref(self.impl_).entity())
        e._entity_pointer = self
        return e

    def __dealloc__(self):

        self.impl_.reset()

    property codim:
        def __get__(self):
            return deref(self.impl_).codimension

    property entity:
        """Return the entity associated with the EntityPointer"""
        def __get__(self):
            return self._entity()


cdef class EntityPointer1:

    cpdef Entity1 _entity(self):
        cpdef Entity1 e = Entity1()
        e.impl_ = address(deref(self.impl_).entity())
        e._entity_pointer = self
        return e

    def __dealloc__(self):

        self.impl_.reset()

    property codim:
        def __get__(self):
            return deref(self.impl_).codimension

    property entity:
        """Return the entity associated with the EntityPointer"""
        def __get__(self):
            return self._entity()

cdef class EntityPointer2:

    cpdef Entity2 _entity(self):
        cpdef Entity2 e = Entity2()
        e.impl_ = address(deref(self.impl_).entity())
        e._entity_pointer = self
        return e

    def __dealloc__(self):

        self.impl_.reset()

    property codim:
        def __get__(self):
            return deref(self.impl_).codimension

    property entity:
        """Return the entity associated with the EntityPointer"""
        def __get__(self):
            return self._entity()

