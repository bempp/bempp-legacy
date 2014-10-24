<%
from sys import version
from bempp_operators import dtypes, compatible_dtypes, bops
ifloop = lambda x: "if" if getattr(x, 'index', x) == 0 else "elif"
division = '__div__' if int(version[0]) < 3 else '__truediv__'
%>\
from bempp.options cimport Options
from bempp.space.space cimport Space

# Declares complex type explicitly.
# Cython 0.20 will fail if templates are nested more than three-deep,
# as in shared_ptr[ c_Space[ complex[float] ] ]
cdef extern from "bempp/space/types.h":
% for ctype in dtypes.values():
%     if 'complex'  in ctype:
    ctypedef struct ${ctype}
%     endif
% endfor

% for name, op in {'multiply': '*', 'divid': '/'}.items():
cdef object _scalar_${name}(BoundaryOperator self,
            BoundaryOperator result, scalar):
    from numpy import iscomplex
    cdef:
        float asfloat
        double asdouble
        float complex ascfloat
        double complex ascdouble
    if iscomplex(scalar):
        if self.result_type == 'complex64':
            ascfloat = scalar
            result.impl_.assign(self.impl_ ${op} ascfloat)
        elif self.result_type == 'complex128':
            ascdouble = scalar
            result.impl_.assign(self.impl_ ${op} ascdouble)
        else:
            raise TypeError(
               "Unable to scale operator using object %s" % str(scalar))
    else:
        if self.result_type in ['complex64', 'float']:
            asfloat = scalar
            result.impl_.assign(self.impl_ ${op} asfloat)
        else:
            asdouble = scalar
            result.impl_.assign(self.impl_ ${op}  asdouble)
    return result
% endfor


cdef class BoundaryOperator:
    """ Holds a reference to a boundary operator """
    def __init__(self, basis_type=None, result_type=None):
        from numpy import dtype

        # Simple way to check basis and result type are compatible
        ops = Options(basis_type=basis_type, result_type=result_type)

% for i, (pybasis, cybasis) in enumerate(dtypes.items()):
%     for j, (pyresult, cyresult) in enumerate(dtypes.items()):
%         if pyresult in compatible_dtypes[pybasis]:
        if ops.basis_type == "${pybasis}" and ops.result_type == "${pyresult}":
            self.impl_.set${cybasis}${cyresult}()
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
            if not self.impl_.valid_${variable}():
                return None
            cdef Space result = Space.__new__(Space)
            result.impl_ = self.impl_.${variable}()
            return result
% endfor

    property label:
        def __get__(self):
            return self.impl_.label()

% for opname, op in {'add': '+', 'sub': '-'}.items():
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
        if not isinstance(self, BoundaryOperator):
            return other.__mul__(self)

        cdef BoundaryOperator res = ${'\\'}
            BoundaryOperator(self.basis_type, self.result_type)
        if isinstance(other, BoundaryOperator):
            arecompatible =                                     ${'\\'}
                (<BoundaryOperator> self).basis_type            ${'\\'}
                    == (<BoundaryOperator> other).basis_type    ${'\\'}
                and (<BoundaryOperator> self).result_type       ${'\\'}
                    == (<BoundaryOperator> other).result_type
            if not arecompatible:
                raise TypeError("Operators are not compatible")
            res.impl_.assign(
                (<BoundaryOperator> self).impl_
                * (<BoundaryOperator> other).impl_
            )
        else:
            _scalar_multiply(self, res, other)
        return res

    def ${division}(BoundaryOperator self, other):
        cdef BoundaryOperator res = BoundaryOperator.__new__(BoundaryOperator)
        if isinstance(other, BoundaryOperator):
            raise TypeError("Cannot divide by an operator")
        _scalar_divid(self, res, other)
        return res
