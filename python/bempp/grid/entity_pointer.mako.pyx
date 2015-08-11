<%
codims = [('0','codim_zero'),('1','codim_one'),('2','codim_two')]
%>

from cython cimport address
from cython.operator cimport dereference as deref

% for (codim,codim_template) in codims:

from bempp.grid.entity cimport Entity${codim}

cdef class EntityPointer${codim}:

    cpdef Entity${codim} _entity(self):
        cpdef Entity${codim} e = Entity${codim}()
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


% endfor     
