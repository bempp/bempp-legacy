<%
from sys import version
from bempp_operators import dtypes, compatible_dtypes, bops
from data_types import scalar_cython_type
%>
from bempp.space.space cimport Space
from discrete_boundary_operator cimport DiscreteBoundaryOperator
from discrete_boundary_operator cimport DiscreteBoundaryOperatorBase
from bempp.utils.byte_conversion import convert_to_bytes
from bempp.assembly.grid_function cimport GridFunction

from numpy cimport dtype

from bempp.utils cimport complex_float,complex_double
import numpy as np
cimport numpy as np


cdef class BoundaryOperatorBase:
    """

    A boundary operator is an abstract representation of an operator
    acting on the boundary of a mesh. It can be multiplied with GridFunction
    objects and provides an operation
    `weak_form` that returns a discrete representation of the boundary
    operator.

    Notes
    -----
    Boundary operators are not constructed directly but through factory
    functions (see for example bempp.operators)

    """
    


    def __cinit__(self):
        pass

    def __init__(self):
        pass

    property basis_type:
        def __get__(self):
            return self._basis_type

    property result_type:
        def __get__(self):
            return self._result_type

    def __add__(self,BoundaryOperatorBase other):

       return _SumBoundaryOperator(self,other)

    def __mul__(self, object x):

        if not isinstance(self,BoundaryOperatorBase):
            return x*self

        if np.isscalar(x):
            return _ScaledBoundaryOperator(self,x)

        if isinstance(x,GridFunction):
            return self._apply_grid_function(x)

        raise ValueError("Multiplication not supported for type"+str(type(x)))

    def __neg__(self):

        return _ScaledBoundaryOperator(self,-1.0)

    def __sub__(self,BoundaryOperatorBase other):

        return self.__add__(-other)

    cpdef DiscreteBoundaryOperatorBase weak_form(self):
        """

        This method returns an assembled representation
        of a boundary operator that allows matrix-vector
        multiplication.

        Returns
        -------
        discrete_operator : bempp.assembly.DiscreteBoundaryOperatorBase
            A discrete representation of the boundary operator.

        """

        raise NotImplementedError("Method not implemented")

    def _apply_grid_function(self,GridFunction g):

        raise NotImplementedError("Method not implemented")

cdef class BoundaryOperator(BoundaryOperatorBase):
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

    def _apply_grid_function(self,GridFunction g):
        if not(self.domain.is_compatible(g.space)):
            raise ValueError("Spaces do not match")
        if not self.result_type==g.result_type:
            raise ValueError("Result types do not match")

        op_w =  self.weak_form()
        coeffs = g.coefficients
        result_projections = (op_w*coeffs)
        return GridFunction(g.parameter_list,self.range,dual_space=self.dual_to_range,
                projections=result_projections)

    cpdef DiscreteBoundaryOperatorBase weak_form(self):
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


    property label:
        def __get__(self):
            return self.impl_.label().decode("UTF-8")

cdef class _ScaledBoundaryOperator(BoundaryOperatorBase):
    cdef BoundaryOperatorBase op
    cdef string _label
% for pyvalue,cyvalue in dtypes.items():
    cdef ${scalar_cython_type(cyvalue)} _alpha_${pyvalue}
% endfor

    def __cinit__(self,BoundaryOperatorBase op,object alpha):

        self._basis_type = op._basis_type
        self._result_type = op._result_type
        self.op = op
% for pyvalue in dtypes:
        if  self._result_type == "${pyvalue}":
            self._alpha_${pyvalue} = alpha
            self._label = bytes(str(alpha)+"*"+op.label,"UTF-8")
% endfor

    def __init__(self,BoundaryOperator op, object alpha):
        pass

    property label:
        def __get__(self):
            return self._label.decode("UTF-8")

    cpdef DiscreteBoundaryOperatorBase weak_form(self):
% for pyvalue in dtypes:
        if self._result_type == "${pyvalue}":
            return self._alpha_${pyvalue}*self.op.weak_form()
% endfor
        raise ValueError("Unknown result type")


cdef class _SumBoundaryOperator(BoundaryOperatorBase):
    cdef BoundaryOperatorBase op1
    cdef BoundaryOperatorBase op2
    cdef string _label

    def __cinit__(self,BoundaryOperatorBase op1, BoundaryOperatorBase op2):

        if not (op1._basis_type==op2._basis_type and op1._result_type==op2._result_type):
            raise ValueError("Operators have different types.")

        self._basis_type = op1._basis_type
        self._result_type = op2._result_type

        self.op1 = op1
        self.op2 = op2

        self._label = bytes("("+op1.label+"+"+op2.label+")","UTF-8")

    def __init__(self,BoundaryOperatorBase op1, BoundaryOperatorBase op2):

        pass

    property label:

        def __get__(self):
            return self._label.decode("UTF-8")

    cpdef DiscreteBoundaryOperatorBase weak_form(self):

        return self.op1.weak_form()+self.op2.weak_form()






