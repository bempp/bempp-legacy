<%
codims = [('0','codim_zero'),('1','codim_one'),('2','codim_two')]
%>

from bempp.utils cimport unique_ptr
from cython.operator cimport dereference as deref
from cython cimport address
from bempp.grid.codim_template cimport codim_zero,codim_one,codim_two

from bempp.grid.entity_iterator cimport c_EntityIterator

% for (codim,codim_template) in codims:
from bempp.grid.entity_iterator cimport EntityIterator${codim}
% endfor

from bempp.utils.eigen cimport eigen_matrix_to_np_float64, eigen_matrix_to_np_int
from bempp.utils cimport Matrix, Vector
from libcpp.vector cimport vector

import numpy as _np
cimport numpy as _np

cdef class GridView:
    """GridView information

       A GridView object contains a view on the entities of a grid.
    """


    def __cinit__(self):
        self._raw_data_is_computed = False

    def __dealloc__(self):
        self.impl_.reset()

    def __init__(self):
        pass

    cpdef size_t entity_count(self,int codim):
        """Return the number of entities of the given codim."""
        return deref(self.impl_).entityCount(codim)

% for (codim, codim_template) in codims:
    cpdef EntityIterator${codim} _entity_iterator${codim}(self):
        """Return an entity iterator for entities of codim ${codim}."""
        cdef:
            EntityIterator${codim} it = EntityIterator${codim}()
            unique_ptr[c_EntityIterator[${codim_template}]] c_it = ${'\\'}
                deref(self.impl_).entityIterator${codim}()
        it.impl_.swap(c_it)
        return it
% endfor

    cdef void _compute_raw_element_data(self):
        if self._raw_data_is_computed: return

        cdef:
            Matrix[double] vertices
            Matrix[int] elements
            Matrix[char] aux_data

        deref(self.impl_).getRawElementData(vertices,elements,aux_data,self._domain_indices)

        self._vertices = eigen_matrix_to_np_float64(vertices)
        self._elements = eigen_matrix_to_np_int(elements)[:-1,:] # Last row not needed for triangular grids

        return


    def entity_iterator(self,codim):
        """Return iterator for entities of given codim."""
        if codim == 0:
            return self._entity_iterator0()
        elif codim == 1:
            return self._entity_iterator1()
        elif codim == 2:
            return self._entity_iterator2()
        else:
            raise RuntimeError()

    property dim:
        """" Dimension of the grid. """
        def __get__(self):
            return deref(self.impl_).dim()

    property dim_world:
        """ Dimension of the space containing the grid. """
        def __get__(self):
            return deref(self.impl_).dimWorld()

    property vertices:
        """ Return the vertices of the grid. """

        def __get__(self):
            self._compute_raw_element_data()
            return self._vertices

    property elements:
        """ Return the elements of the grid. """

        def __get__(self):
            self._compute_raw_element_data()
            return self._elements

    property domain_indices:
        """ Return the domain indices of the elements. """

        def __get__(self):
            self._compute_raw_element_data()
            return self._domain_indices



cdef GridView _grid_view_from_unique_ptr(unique_ptr[c_GridView]& c_view):
    cdef GridView view = GridView()
    view.impl_.swap(c_view)
    return view
