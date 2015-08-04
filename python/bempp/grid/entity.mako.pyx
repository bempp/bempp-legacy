<%
codims = [('0','codim_zero'),('1','codim_one'),('2','codim_two')]
%>

from cython.operator cimport dereference as deref
from cython cimport address
from bempp.utils cimport unique_ptr
from bempp.grid.codim_template cimport codim_zero,codim_one,codim_two
from bempp.grid.entity_pointer cimport c_EntityPointer, EntityPointer0

% for (codim,codim_template) in codims:

from bempp.grid.geometry cimport Geometry${codim}

cdef extern from "bempp/grid/py_entity_helper.hpp" namespace "Bempp":
    cdef unique_ptr[c_EntityPointer[codim_zero]] py_get_father_from_entity(const c_Entity[codim_zero]&)

cdef class Entity${codim}:

    cpdef size_t level(self):
        return deref(self.impl_).level()

    cpdef Geometry${codim} _geometry(self):
        cdef Geometry${codim} g = Geometry${codim}()
        g.impl_= address(deref(self.impl_).geometry())
        g._entity = self
        return g

    property geometry:
        """Return the geometry of the entity"""
        def __get__(self):
            return self._geometry()

%if codim=='0':
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

%endif


% endfor


