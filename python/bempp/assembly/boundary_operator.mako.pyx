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

    property label:
        def __get__(self):
            return self.impl_.label()

% for opname, op in {'add': '+', 'sub': '-'}.iteritems():
    def __${opname}__(self, other):
        if not (
            isinstance(self, BoundaryOperator)
            and isinstance(other, BoundaryOperator)
        ):
            raise TypeError("Incorrect types in ${opname}")
        cdef BoundaryOperator res = BoundaryOperator.__new__(BoundaryOperator)
        res.impl_.assign(
            (<BoundaryOperator> self).impl_
            ${op} (<BoundaryOperator> other).impl_
        )
        return res
% endfor

    def __mul__(self, other):
        cdef BoundaryOperator res = BoundaryOperator.__new__(BoundaryOperator)
        cdef double asdouble
        if (
            isinstance(self, BoundaryOperator)
            and isinstance(other, BoundaryOperator)
        ):
            res.impl_.assign(
                (<BoundaryOperator> self).impl_
                * (<BoundaryOperator> other).impl_
            )
        elif not isinstance(self, BoundaryOperator):
            return other.__mul__(self)
        elif isinstance(other, complex):
            res.impl_.assign(
                    (<BoundaryOperator> self).impl_ * (<complex>other))
        else:
            try:
                asdouble = other
                res.impl_.assign((<BoundaryOperator> self).impl_ * asdouble)
            except:
                raise TypeError("Incorrect types in multiplication")
        return res

    def __div__(BoundaryOperator self, other):
        cdef BoundaryOperator res = BoundaryOperator.__new__(BoundaryOperator)
        cdef double asdouble
        if isinstance(other, BoundaryOperator):
            raise TypeError("Cannot divide by an operator")
        elif isinstance(other, complex):
            res.impl_.assign(
                    (<BoundaryOperator> self).impl_ * (<complex>other))
        else:
            try:
                asdouble = other
                res.impl_.assign((<BoundaryOperator> self).impl_ * asdouble)
            except:
                raise TypeError("Incorrect types in division")
        return res
