<%
codims = [('0','codim_zero'),('1','codim_one'),('2','codim_two')]
%>

from bempp.utils cimport unique_ptr
from bempp.utils cimport catch_exception
from cython.operator cimport dereference as deref
from cython cimport address
from libcpp cimport bool as cbool


from bempp.grid.entity cimport c_Entity
from bempp.grid.entity_pointer cimport c_EntityPointer

% for (codim,codim_template) in codims:

from bempp.grid.entity cimport Entity${codim}
from bempp.grid.entity_pointer cimport EntityPointer${codim}

cdef class EntityIterator${codim}:

    cdef cbool finished(self):
        return deref(self.impl_).finished()

    cdef void c_next(self) except + :
        deref(self.impl_).next()

    cdef EntityPointer${codim} _frozen(self):
        cdef EntityPointer${codim} ep = EntityPointer${codim}()
        cdef unique_ptr[c_EntityPointer[${codim_template}]] c_ep = deref(self.impl_).frozen()
        ep.impl_.swap(c_ep)
        return ep

    def __next__(self):
        cdef EntityPointer${codim} ep
        if not self.finished():
            ep = self._frozen()
            self.c_next()
            return ep.entity
        else:
            raise StopIteration()

    def __iter__(self):
        return self

% endfor





       




        


