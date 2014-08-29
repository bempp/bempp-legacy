<%
from space import dtypes, spaces, flatten

def declare_class(text):
    return 'c_{0} "Bempp::{0}"[BASIS]'.format(text)
%>
from libcpp cimport complex as ccomplex
from bempp.utils cimport shared_ptr
from bempp.grid.grid cimport Grid, c_Grid

cdef extern from "bempp/space/space.hpp":
    cdef cppclass c_Space "Bempp::Space" [BASIS]:
        c_Space(const shared_ptr[c_Grid]&)
        c_Space(const c_Space[BASIS]&)

# Declares complex type explicitly.
# Cython 0.20 will fail if templates are nested more than three-deep,
# as in shared_ptr[ c_Space[ complex[float] ] ]
cdef extern from "bempp/space/types.h":
% for ctype in dtypes.itervalues():
%     if 'complex'  in ctype:
    ctypedef struct ${ctype}
%     endif
% endfor



% for space, description in flatten(spaces):
cdef extern from "bempp/space/${space.header}.hpp":
    cdef cppclass ${space.class_name | declare_class}:
%   if description['implementation'] == 'grid_only':
        ${'c_' + space.class_name}(const shared_ptr[c_Grid] &_grid)
%   endif
% endfor

cdef class Space:
    cdef:
        readonly Grid grid
        # For simplicity, we define shared pointers to all possible space types
        # There are not that many, so this should not be a problem
% for pytype, cytype in dtypes.iteritems():
        shared_ptr[c_Space[${cytype}]] impl_${pytype}
% endfor

# Now we define derived types for each space.
# This is a flat hierarchy. It does not attempt to redeclare the C++ hierarchy.
% for space, description in flatten(spaces):
cdef class ${space.class_name}(Space):
    pass
% endfor
