<%
from sys import version
from bempp_operators import dtypes, compatible_dtypes, bops
from data_types import scalar_cython_type
%>
from bempp.space.space cimport Space
from discrete_boundary_operator cimport DiscreteBoundaryOperator
from discrete_boundary_operator cimport DenseDiscreteBoundaryOperator
from discrete_boundary_operator cimport DiscreteBoundaryOperatorBase
from discrete_boundary_operator import SparseDiscreteBoundaryOperator
from discrete_boundary_operator import ZeroDiscreteBoundaryOperator
from bempp.utils.byte_conversion import convert_to_bytes
from bempp.assembly.grid_function cimport GridFunction
from bempp.utils cimport shared_ptr,static_pointer_cast
from bempp.utils import combined_type
from bempp.space.space cimport Space

from numpy cimport dtype

from bempp.utils cimport complex_float,complex_double
import numpy as np
cimport numpy as np

np.import_array()

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

    property domain:

        def __get__(self):
            raise ValueError("Method not implemented.")

    property range:

        def __get__(self):
            raise ValueError("Method not implemented.")

    property dual_to_range:

        def __get__(self):
            raise ValueError("Method not implemented.")

    def __add__(self,BoundaryOperatorBase other):

       return _SumBoundaryOperator(self,other)

    def __mul__(self, object x):

        if not isinstance(self,BoundaryOperatorBase):
            return x*self

        if np.isscalar(x):
            return _ScaledBoundaryOperator(self,x)

        if isinstance(x,GridFunction):
            return self._apply_grid_function(x)

        raise NotImplementedError("Multiplication not supported for type"+str(type(x)))

    def __neg__(self):

        return _ScaledBoundaryOperator(self,-1.0)

    def __sub__(self,BoundaryOperatorBase other):

        return self.__add__(-other)

    def weak_form(self):
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

cdef class GeneralBoundaryOperator(BoundaryOperatorBase):
    """ Holds a reference to a boundary operator """


    def __cinit__(self, basis_type=None, result_type=None):
        pass


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

        op_w =  self.weak_form()
        coeffs = g.coefficients
        result_projections = (op_w*coeffs)
        return GridFunction(self.range,dual_space=self.dual_to_range,
                projections=result_projections,parameter_list=g.parameter_list)

    def weak_form(self):
        cdef DiscreteBoundaryOperator dbop = DiscreteBoundaryOperator()
        
% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:

        if self.basis_type=="${pybasis}" and self.result_type=="${pyresult}":
            dbop._impl_${pyresult}_.assign(_boundary_operator_variant_weak_form[${cybasis},${cyresult}](self.impl_))
            dbop._dtype = self.result_type
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
    cdef BoundaryOperatorBase _op
    cdef string _label
    cdef object _alpha

    def __cinit__(self,BoundaryOperatorBase op,object alpha):
        pass

    def __init__(self,BoundaryOperatorBase op, object alpha):
        
        self._alpha = 1.0*alpha # make sure it is not integer
        self._basis_type = op._basis_type
        self._result_type = combined_type(np.dtype(type(self._alpha)),op.result_type)
        self._op = op

        self._label = (str(alpha)+"*"+self._op.label).encode("UTF-8")

    property label:
        def __get__(self):
            return self._label.decode("UTF-8")

    def weak_form(self):
        return self._alpha*self._op.weak_form()

    def _apply_grid_function(self,GridFunction g):

        return self._alpha*(self._op*g)


    property domain:

        def __get__(self):
            return self._op.domain

    property range:

        def __get__(self):
            return self._op.range

    property dual_to_range:

        def __get__(self):
            return self._op.dual_to_range
    

cdef class _SumBoundaryOperator(BoundaryOperatorBase):
    cdef BoundaryOperatorBase _op1
    cdef BoundaryOperatorBase _op2
    cdef string _label

    def __cinit__(self,BoundaryOperatorBase op1, BoundaryOperatorBase op2):
        pass


    def __init__(self,BoundaryOperatorBase op1, BoundaryOperatorBase op2):

        self._basis_type = combined_type(op1._basis_type,op2._basis_type)
        self._result_type = combined_type(op1._result_type,op2._result_type)

        if not (op1.range.is_compatible(op2.range) and 
                op1.dual_to_range.is_compatible(op2.dual_to_range) and
                op1.domain.is_compatible(op2.domain)):
            raise ValueError("Incompatible Spaces")

        self._op1 = op1
        self._op2 = op2

        self._label = ("("+op1.label+"+"+op2.label+")").encode("UTF-8")

    property label:

        def __get__(self):
            return self._label.decode("UTF-8")

    def weak_form(self):

        return self._op1.weak_form()+self._op2.weak_form()

    def _apply_grid_function(self, GridFunction g):

        return self._op1*g+self._op2*g

    property domain:

        def __get__(self):
            return self._op1.domain

    property range:

        def __get__(self):
            return self._op1.range

    property dual_to_range:

        def __get__(self):
            return self._op1.dual_to_range

cdef class SparseBoundaryOperator(GeneralBoundaryOperator):
    """ Holds a reference to a boundary operator """



    def __cinit__(self, basis_type=None, result_type=None):
        pass


    def __init__(self, basis_type=None, result_type=None):

        super(SparseBoundaryOperator,self).__init__(basis_type,result_type)

    def weak_form(self):

        import scipy.sparse

        cdef DiscreteBoundaryOperator discrete_operator = super(SparseBoundaryOperator,self).weak_form()
        res = None
% for pyvalue,cyvalue in dtypes.items():

        if discrete_operator.dtype=="${pyvalue}":
            res = py_get_sparse_from_discrete_operator[${cyvalue}](discrete_operator._impl_${pyvalue}_)
% endfor

        if res is None:
            raise ValueError("Unknown data type")

        return SparseDiscreteBoundaryOperator(res)


    property domain:

        def __get__(self):
            return self._domain

    property range:

        def __get__(self):
            return self._range

    property dual_to_range:

        def __get__(self):
            return self._dual_to_range

cdef class DenseBoundaryOperator(GeneralBoundaryOperator):


    def __cinit__(self, GeneralBoundaryOperator op):
        pass


    def __init__(self,GeneralBoundaryOperator op):

        super(DenseBoundaryOperator,self).__init__(basis_type=op.basis_type,result_type=op.result_type)
        self.impl_ = op.impl_

    def weak_form(self):
        cdef DenseDiscreteBoundaryOperator dbop = DenseDiscreteBoundaryOperator()
        
% for pybasis,cybasis in dtypes.items():
%     for pyresult,cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:

        if self.basis_type=="${pybasis}" and self.result_type=="${pyresult}":
            dbop._impl_${pyresult}_.assign(_boundary_operator_variant_weak_form[${cybasis},${cyresult}](self.impl_))
            dbop._dtype = self.result_type
            return dbop
%          endif
%      endfor
% endfor
        raise ValueError("Incompatible basis and result types") 

cdef class ZeroBoundaryOperator(BoundaryOperatorBase): 

    cdef Space _domain
    cdef Space _range
    cdef Space _dual_to_range

    def __cinit__(self,Space domain, Space range, Space dual_to_range):
        pass

    def __init__(self,Space domain, Space range, Space dual_to_range):

        self._domain = domain
        self._range = range
        self._dual_to_range = dual_to_range
        self._basis_type = np.dtype('float64')
        self._result_type = np.dtype('float64')

    property basis_type:
        def __get__(self):
            return self._basis_type

    property result_type:
        def __get__(self):
            return self._result_type

    property domain:

        def __get__(self):
            return self._domain

    property range:

        def __get__(self):
            return self._range

    property dual_to_range:

        def __get__(self):
            return self._dual_to_range

    def weak_form(self):
        """

        This method returns an assembled representation
        of a boundary operator that allows matrix-vector
        multiplication.

        Returns
        -------
        discrete_operator : bempp.assembly.DiscreteBoundaryOperatorBase
            A discrete representation of the boundary operator.

        """

        return ZeroDiscreteBoundaryOperator(self._domain.global_dof_count,
                self._dual_to_range.global_dof_count)

    def _apply_grid_function(self,GridFunction g):
        return GridFunction(self._domain,
                coefficients=np.zeros(self._domain.global_dof_count,
                    dtype='float64'))

cdef class BlockedBoundaryOperator:

    cdef object _operators
    cdef object _domain_spaces
    cdef object _range_spaces 
    cdef object _dual_to_range_spaces

    def __cinit__(self,domain_spaces,range_spaces,dual_to_range_spaces):
        pass

    def __init__(self,domain_spaces,range_spaces,dual_to_range_spaces):

        if not len(dual_to_range_spaces)==len(range_spaces):
            raise ValueError("Dimension of 'dual_to_range_spaces' must be identical to dimension of 'range_spaces'")

        self._domain_spaces = domain_spaces
        self._range_spaces = range_spaces
        self._dual_to_range_spaces = dual_to_range_spaces

        self._operators = np.empty(self.ndims,dtype=np.object)

        for i,range_space in enumerate(range_spaces):
            for j,domain_space in enumerate(domain_spaces):
                self._operators[i,j] = ZeroBoundaryOperator(
                        domain_space,range_space,
                        self.dual_to_range[i])

    property ndims:

        def __get__(self):

            return (len(self.range),len(self.domain))

    property domain:

        def __get__(self):
            return self._domain_spaces

    property range:

        def __get__(self):
            return self._range_spaces
    
    property dual_to_range:

        def __get__(self):
            return self._dual_to_range_spaces

    property shape:

        def __get__(self):
            return self._operators.shape

    def __getitem__(self,key):
        return self._operators[key]

    def __setitem__(self,key,BoundaryOperatorBase op):

        if not self[key].domain.is_compatible(op.domain):
            raise ValueError("Domain space is not compatible")

        if not self[key].range.is_compatible(op.range):
            raise ValueError("Range Space is not compatible")

        if not self[key].dual_to_range.is_compatible(op.dual_to_range):
            raise ValueError("Dual space is not compatible")
        
        self._operators[key] = op

    def __add__(self,BlockedBoundaryOperator other):

        if not self.shape==other.shape:
            raise ValueError("Dimensions not compatible")
        
        sum_operator = BlockedBoundaryOperator(
                self.domain,self.range,self.dual_to_range)

        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                sum_operator[i,j] = self[i,j]+other[i,j]
        return sum_operator
                
    def __mul__(self, object x):

        if not isinstance(self,BlockedBoundaryOperator):
            return x*self

        if np.isscalar(x):
            scaled_blocked_operator = BlockedBoundaryOperator(
                    self.domain,self.range,self._dual_to_range)
            for i in range(self.shape[0]):
                for j in range(self.shape[1]):
                    scaled_blocked_operator[i,j] = x*self[i,j]
            return scaled_blocked_operator

        for g in x:
            if not isinstance(g,GridFunction):
                raise NotImplementedError("Multiplication not supported for type"+str(type(x)))

        res = []
        for row in self.shape[0]:
            res_fun = GridFunction(self.range[0],
                    coefficients=np.zeros(self.range[0].global_dof_count))
            for col in self.shape[1]:
                res_fun = res_fun+self[row,col]*x[col]
            res.append(res_fun)
        return res

    def __neg__(self):

        return self.__mul__(-1.0)

    def __sub__(self,BlockedBoundaryOperator other):

        return self.__add__(-other)

    def weak_form(self):
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


