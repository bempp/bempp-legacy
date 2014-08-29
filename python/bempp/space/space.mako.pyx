<%
from space import dtypes, spaces, flatten


def long_names(spaces, current_name="space"):
    for key, values in spaces.iteritems():
        long_name = current_name + "." + key
        iterator = values if key is None else long_names(values, long_name)
        for value in iterator:
            yield long_name, value
%>


cdef class Space:
    """ Space of functions defined on a grid

        The exact space depends on the input.
    """
    def __init__(self, Grid grid not None):
        self.grid = grid

    property dtype:
        """ Precision and kind of this space """
        def __get__(self):
            from numpy import dtype
% for pyname, cython in dtypes.iteritems():
            if self.impl_${pyname}.get() is not NULL:
                return dtype('${pyname}')
% endfor
            return None

% for space, description in flatten(spaces):
cdef class ${space.class_name}(Space):
    """ ${description['doc']}

        Parameters:
        ----------

        grid : Grid
            Grid over which to discrtize the space

        dtype : numpy.dtype
            Type of the functions acting on the space
    """
    def __init__(self, Grid grid not None, dtype, *args, **kwargs):
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
                new ${'c_' + space.class_name}[${cytype}](grid.impl_)
            )
%       endif
%    endfor

% endfor
