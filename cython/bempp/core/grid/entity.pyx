from cython.operator cimport dereference as deref
from cython cimport address
from bempp.core.utils cimport unique_ptr
from bempp.core.grid.codim_template cimport codim_zero,codim_one,codim_two
from bempp.core.grid.entity_pointer cimport c_EntityPointer, EntityPointer0
from bempp.core.grid.entity_iterator cimport EntityIterator1, EntityIterator2
from bempp.core.grid.entity_iterator cimport c_EntityIterator
from bempp.core.grid.geometry cimport Geometry0
from bempp.core.grid.geometry cimport Geometry1
from bempp.core.grid.geometry cimport Geometry2

cdef extern from "bempp/core/grid/py_entity_helper.hpp" namespace "Bempp":
    cdef unique_ptr[c_EntityPointer[codim_zero]] py_get_father_from_entity(const c_Entity[codim_zero]&)
    cdef unique_ptr[c_EntityIterator[ITERATOR_CODIM]] py_get_entity_iterator_from_entity[ITERATOR_CODIM, ENTITY_CODIM]\
        (const c_Entity[ENTITY_CODIM])

cdef class Entity0:

    cpdef size_t level(self):
        return deref(self.impl_).level()

    cpdef Geometry0 _geometry(self):
        cdef Geometry0 g = Geometry0()
        g.impl_= address(deref(self.impl_).geometry())
        g._entity = self
        return g

    def sub_entity_iterator(self, int codim):

        cdef EntityIterator1 it1 = EntityIterator1()
        cdef EntityIterator2 it2 = EntityIterator2()

        cdef unique_ptr[c_EntityIterator[codim_one]] ptr1
        cdef unique_ptr[c_EntityIterator[codim_two]] ptr2


        if codim == 1:
            ptr1 = py_get_entity_iterator_from_entity[codim_one, codim_zero](deref(self.impl_))
            it1.impl_.swap(ptr1)
            return it1
        elif codim == 2:
            ptr2 = py_get_entity_iterator_from_entity[codim_two, codim_zero](deref(self.impl_))
            it2.impl_.swap(ptr2)
            return it2
        else:
            raise ValueError("Wrong codimension.")



    property geometry:
        """Return the geometry of the entity"""
        def __get__(self):
            return self._geometry()

    property domain:
        """Return the domain index of the entity"""
        def __get__(self):
            return deref(self.impl_).domain()

    property is_leaf:
        """Check if an element is a leaf"""
        def __get__(self):
            return deref(self.impl_).isLeaf()

    property has_father:
        """Check if entity has a father in the grid"""
        def __get__(self):
            return deref(self.impl_).hasFather()

    property father:
        """Return father entity"""
        def __get__(self):
            if not self.has_father:
                raise ValueError("Entity has no father.")
            cdef EntityPointer0 ep = EntityPointer0()
            cdef unique_ptr[c_EntityPointer[codim_zero]] c_ep = py_get_father_from_entity(deref(self.impl_))
            ep.impl_.swap(c_ep)
            return ep.entity


cdef class Entity1:

    cpdef size_t level(self):
        return deref(self.impl_).level()

    cpdef Geometry1 _geometry(self):
        cdef Geometry1 g = Geometry1()
        g.impl_= address(deref(self.impl_).geometry())
        g._entity = self
        return g

    property geometry:
        """Return the geometry of the entity"""
        def __get__(self):
            return self._geometry()


cdef class Entity2:

    cpdef size_t level(self):
        return deref(self.impl_).level()

    cpdef Geometry2 _geometry(self):
        cdef Geometry2 g = Geometry2()
        g.impl_= address(deref(self.impl_).geometry())
        g._entity = self
        return g

    property geometry:
        """Return the geometry of the entity"""
        def __get__(self):
            return self._geometry()

