""" Wrappers for all types of C++ spaces """
<%
from data_types import dtypes
from space import spaces
ifloop = lambda x: 'if' if loop.index == 0 else 'elif'
%>
from bempp.grid.grid cimport Grid
from cython.operator cimport dereference as deref
from libcpp cimport bool as cbool


cdef class Space:
    """ Space of functions defined on a grid

        The exact space depends on the input.
    """
    def __init__(self, Grid grid not None):
        super(Space, self).__init__()

    property dtype:
        """ Precision and kind of this space """
        def __get__(self):
            from numpy import dtype
            return dtype(self.impl_.dtype());

    property codomainDimension:
        """Number of components of values of functions in this space (e.g. 1 for scalar functions)"""
        def __get__(self):
            return self.impl_.codomainDimension()

    property grid:
        def __get__(self):
            cdef Grid result = Grid.__new__(Grid)
            result.impl_ = self.impl_.grid()
            return result

    cpdef cbool is_compatible(self,Space other):
        return self.impl_.isCompatible(other.impl_)

% for class_name, description in spaces.items():
cdef class ${class_name}(Space):
    """ ${description['doc']}

        Parameters:
        ----------

        grid : Grid
            Grid over which to discretize the space

        dtype : numpy.dtype
            Type of the functions acting on the space

% if description['implementation'] == 'polynomial':
        order : int
            Order of the polynomial. 
% endif
    """
    def __init__(self, Grid grid not None, dtype,
% if description['implementation'] == 'polynomial':
                 order,
% endif
                 **kwargs):
        from numpy import dtype as np_dtype
        super(${class_name}, self).__init__(grid)

        dtype = np_dtype(dtype)
        if dtype not in ${list(dtypes.keys())}:
                raise TypeError("Unexpected basis type")
%    for pytype, cytype in dtypes.items():
        if dtype == "${pytype}":
%       if description['implementation'] == 'grid_only':
            self.impl_.set( shared_ptr[c_Space[${cytype}]](
                <c_Space[${cytype}]*>
                new ${'c_' + class_name}[${cytype}](grid.impl_)
            ))
%       elif description['implementation'] == 'polynomial':
            self.order = order
            self.impl_.set( shared_ptr[c_Space[${cytype}]](
                <c_Space[${cytype}]*> new ${'c_' + class_name}[${cytype}](
                    grid.impl_, <int> self.order
                )
            ))
%       endif
%    endfor
% endfor
