""" Wrappers for all types of C++ spaces """
<%
from data_types import dtypes
from space import spaces
ifloop = lambda x: 'if' if loop.index == 0 else 'elif'
%>
from bempp.grid.grid cimport Grid
from cython.operator cimport dereference as deref


cdef class Space:
    """ Space of functions defined on a grid

        The exact space depends on the input.
    """
    def __init__(self, Grid grid not None):
        super(Space, self).__init__()

    def is_compatible(self, Space other):
        """ True if other is compatible with this space

            Two spaces are compatible if:

            * They have the same dtype
            * Their global degress of freedom agree
        """
        if other is None:
            return False
        return self.impl_.isCompatible(other.impl_)

    property dtype:
        """ Precision and kind of this space """
        def __get__(self):
            from numpy import dtype
            return dtype(self.impl_.dtype());

    property grid:
        def __get__(self):
            cdef Grid result = Grid.__new__(Grid)
            result.impl_ = self.impl_.grid()
            return result

    def __richcmp__(Space self, Space other not None, int op):
        if op != 2:
            raise AttributeError("Incorrect operator")
        return self.impl_.is_same(other.impl_)


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
            Order of the polynomial. Defaults to 2.
% endif
    """
    def __init__(self, Grid grid not None, dtype,
% if description['implementation'] == 'polynomial':
                 order=2,
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
