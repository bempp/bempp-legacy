""" Wrappers for all types of C++ spaces """
<%
from space import dtypes, spaces
ifloop = lambda x: 'if' if loop.index == 0 else 'elif'
%>
from bempp.grid.grid cimport Grid
from cython.operator cimport dereference as deref


cdef class Space:
    """ Space of functions defined on a grid

        The exact space depends on the input.
    """
    def __init__(self, Grid grid not None):
        pass

    property dtype:
        """ Precision and kind of this space """
        def __get__(self):
            from numpy import dtype
% for pyname, cython in dtypes.iteritems():
            ${ifloop(loop)} self.impl_${pyname}.get() is not NULL:
                return dtype('${pyname}')
% endfor
            return None

    property grid:
        """ Grid over which the space is defined """
        def __get__(self):
            if self.dtype is None:
                return None
            cdef Grid grid = Grid.__new__(Grid)
% for pyname in dtypes.iterkeys():
            ${ifloop(loop)} self.dtype == "${pyname}":
                grid.impl_ = deref(self.impl_${pyname}).grid()
% endfor
            return grid

% for class_name, description in spaces.iteritems():
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
            Order of the polynomial. Defaults to 2.
% endif
    """
    def __init__(self, Grid grid not None, dtype,
% if description['implementation'] == 'polynomial':
                 order=2,
% endif
                 **kwargs):
        from numpy import dtype as np_dtype
        super(Space, self).__init__(grid)

        dtype = np_dtype(dtype)
        if dtype not in ${dtypes.keys()}:
                raise TypeError("Unexpected basis type")
%    for pytype, cytype in dtypes.iteritems():
        if dtype == "${pytype}":
%       if description['implementation'] == 'grid_only':
            self.impl_${pytype}.reset(
                <c_Space[${cytype}]*>
                new ${'c_' + class_name}[${cytype}](grid.impl_)
            )
%       elif description['implementation'] == 'polynomial':
            self.order = order
            self.impl_${pytype}.reset(
                <c_Space[${cytype}]*> new ${'c_' + class_name}[${cytype}](
                    grid.impl_, <int> self.order
                )
            )
%       endif
%    endfor
% endfor
