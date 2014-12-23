<%
from sys import version
from bempp_operators import dtypes, compatible_dtypes, bops
ifloop = lambda x: "if" if getattr(x, 'index', x) == 0 else "elif"
division = '__div__' if int(version[0]) < 3 else '__truediv__'
%>\
from bempp.space.space cimport Space
from discrete_boundary_operator cimport DiscreteBoundaryOperator
from numpy cimport dtype

from bempp.utils cimport complex_float,complex_double

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

        self._basis_type = basis_type
        self._result_type = result_type

% for i, (pybasis, cybasis) in enumerate(dtypes.items()):
%     for j, (pyresult, cyresult) in enumerate(dtypes.items()):
%         if pyresult in compatible_dtypes[pybasis]:
        if basis_type == "${pybasis}" and result_type == "${pyresult}":
            self.impl_.set${cybasis}${cyresult}()
%         endif
%     endfor
% endfor

    cpdef DiscreteBoundaryOperator weakForm(self):
        cdef DiscreteBoundaryOperator dbop = DiscreteBoundaryOperator()
        
% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:

        if self.basis_type=="${pybasis}" and self.result_type=="${pyresult}":
            dbop._impl_${pyresult}_.assign(_boundary_operator_variant_weak_form[${cybasis},${cyresult}](self.impl_))
            dbop._value_type = self.result_type
            return dbop
%          endif
%      endfor
% endfor
        raise ValueError("Incompatible basis and result types") 
        

% for variable in ['domain', 'range', 'dual_to_range']:
    property ${variable}:
        def __get__(self):
            if not self.impl_.valid_${variable}():
                return None
            cdef Space result = Space.__new__(Space)
            result.impl_ = self.impl_.${variable}()
            return result
% endfor

    property basis_type:
        def __get__(self):
            return self._basis_type

    property result_type:
        def __get__(self):
            return self._result_type

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
