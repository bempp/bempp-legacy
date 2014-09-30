<%
from bempp_operators import dtypes, compatible_dtypes, bops
ifloop = lambda x: "if" if getattr(x, 'index', x) == 0 else "elif"
%>\
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from bempp.options cimport Options
from bempp.space.space cimport Space

# Declares complex type explicitly.
# Cython 0.20 will fail if templates are nested more than three-deep,
# as in shared_ptr[ c_Space[ complex[float] ] ]
cdef extern from "bempp/space/types.h":
% for ctype in dtypes.itervalues():
%     if 'complex'  in ctype:
    ctypedef struct ${ctype}
%     endif
% endfor

cdef class BoundaryOperator:
    """ Holds a reference to a boundary operator """

    def __init__(self, basis_type=None, result_type=None):
        from numpy import dtype

        # Simple way to check basis and result type are compatible
        ops = Options(basis_type=basis_type, result_type=result_type)

% for i, (pybasis, cybasis) in enumerate(dtypes.iteritems()):
%     for j, (pyresult, cyresult) in enumerate(dtypes.iteritems()):
%         if pyresult in compatible_dtypes[pybasis]:
        if ops.basis_type == "${pybasis}" and ops.result_type == "${pyresult}":
            self.impl_.set[${cybasis}, ${cyresult}]()
%         endif
%     endfor
% endfor

% for variable in ['basis', 'result']:
    property ${variable}_type:
        def __get__(self):
            from numpy import dtype
            return dtype(self.impl_.${variable + "Type"}())
% endfor

% for variable in ['domain', 'range', 'dual_to_range']:
    property ${variable}:
        def __get__(self):
            cdef Space result = Space.__new__(Space)
            try:
                result.impl_ = self.impl_.${variable}()
            except RuntimeError:
                # Unitialized operator
                return None
            else:
                return result
% endfor
