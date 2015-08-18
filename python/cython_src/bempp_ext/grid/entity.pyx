from cython.operator cimport dereference as deref
from cython cimport address
from bempp_ext.utils cimport unique_ptr
from bempp_ext.grid.codim_template cimport codim_zero,codim_one,codim_two
from bempp_ext.grid.entity_pointer cimport c_EntityPointer, EntityPointer0
from bempp_ext.grid.geometry cimport Geometry0
from bempp_ext.grid.geometry cimport Geometry1
from bempp_ext.grid.geometry cimport Geometry2

cdef extern from "bempp_ext/grid/py_entity_helper.hpp" namespace "Bempp":
    cdef unique_ptr[c_EntityPointer[codim_zero]] py_get_father_from_entity(const c_Entity[codim_zero]&)

cdef class Entity0:

    cpdef size_t level(self):
        return deref(self.impl_).level()

    cpdef Geometry0 _geometry(self):
        cdef Geometry0 g = Geometry0()
        g.impl_= address(deref(self.impl_).geometry())
        g._entity = self
        return g

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

