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
from bempp.utils cimport ParameterList

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


    def __cinit__(self, object basis_type, object result_type,
            ParameterList parameters, object is_sparse=False):
        
        from numpy import dtype
        self._basis_type = dtype(basis_type)
        self._result_type = dtype(result_type)
        self._parameters = parameters
        self._is_sparse = is_sparse

    def __init__(self, object basis_type, object result_type,
            ParameterList parameters, object is_sparse=False):
        pass

    def _apply_grid_function(self,GridFunction g):
        if not(self.domain.is_compatible(g.space)):
            raise ValueError("Spaces do not match")

        op_w =  self.weak_form()
        coeffs = g.coefficients
        result_projections = (op_w*coeffs)
        return GridFunction(self.range,dual_space=self.dual_to_range,
                projections=result_projections,parameter_list=g.parameter_list)

    def weak_form(self):

        if self.is_sparse:
            return self._sparse_weak_form()
        elif self.parameter_list.assembly.boundary_assembly_type == 'dense':
            return self._dense_weak_form()
        elif self.parameter_list.assembly.boundary_assembly_type =='hmat':
            return self._hmat_weak_form()
        else:
            return self._default_weak_form()


    def _sparse_weak_form(self):

        # scipy.sparse is needed from C. Imported already here
        # to throw exception if not available.
        import scipy.sparse 

        cdef DiscreteBoundaryOperator discrete_operator = self._default_weak_form()
        res = None
% for pyvalue,cyvalue in dtypes.items():

        if discrete_operator.dtype=="${pyvalue}":
            res = py_get_sparse_from_discrete_operator[${cyvalue}](discrete_operator._impl_${pyvalue}_)
% endfor

        if res is None:
            raise ValueError("Unknown data type")

        return SparseDiscreteBoundaryOperator(res)

    def _dense_weak_form(self):

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

    def _hmat_weak_form(self):

        return self._default_weak_form()

    def _default_weak_form(self):

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

    property parameter_list:

        def __get__(self):
            return self._parameters

    property is_sparse:

        def __get__(self):
            return self._is_sparse

cdef class _ScaledBoundaryOperator(BoundaryOperatorBase):
    cdef BoundaryOperatorBase _op
    cdef object _alpha

    def __cinit__(self,BoundaryOperatorBase op,object alpha):
        pass

    def __init__(self,BoundaryOperatorBase op, object alpha):
        
        self._alpha = 1.0*alpha # make sure it is not integer
        self._basis_type = op._basis_type
        self._result_type = combined_type(np.dtype(type(self._alpha)),op.result_type)
        self._op = op



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

    cdef np.ndarray _operators
    cdef object _ndims
    
    cdef object _range_spaces
    cdef object _domain_spaces
    cdef object _dual_to_range_spaces

    def __cinit__(self,M,N):
        pass

    def __init__(self,M,N):

        self._ndims = (M,N)
        self._operators = np.empty((M,N),dtype=np.object)

        self._domain_spaces = np.empty(N,dtype=np.object)
        self._range_spaces = np.empty(M,dtype=np.object)
        self._dual_to_range_spaces = np.empty(M,dtype=np.object)

    property ndims:

        def __get__(self):
            return self._ndims

    property range_spaces:

        def __get__(self):
            return self._range_spaces

    property domain_spaces:

        def __get__(self):
            return self._domain_spaces

    property dual_to_range_spaces:

        def __get__(self):
            return self._dual_to_range_spaces

    property result_type:

        def __get__(self):
            dt = np.dtype('float64')
            for op in self._operators.ravel():
                if op is not None:
                    dt = combined_type(dt,op.result_type)
            return dt

    property basis_type:

        def __get__(self):
            dt = np.dtype('float64')
            for op in self._operators.ravel():
                if op is not None:
                    dt = combined_type(dt,op.basis_type)
            return dt


    def fill_complete(self):

        for i in range(self.ndims[0]):
            for j in range(self.ndims[1]):
                if (self._domain_spaces[j] is None):
                    raise ValueError("No operator found in column {0}".format(j))
                if (self._range_spaces[i] is None):
                    raise ValueError("No operator found in row {0}".format(i))

                if self._operators[i,j] is None:
                    self._operators[i,j] = ZeroBoundaryOperator(
                            self._domain_spaces[j],self.range_spaces[i],
                            self._dual_to_range_spaces[i])

    def __getitem__(self,key):

        return self._operators[key]

    def __setitem__(self,key,BoundaryOperatorBase op):
        
        i,j = key
        if (self._range_spaces[i] is not None) and (not self._range_spaces[i].is_compatible(op.range)):
                raise ValueError("Range space not compatible")
        if (self._domain_spaces[j] is not None) and (not self._domain_spaces[j].is_compatible(op.domain)):
                raise ValueError("Domain space not compatible")
        if (self._dual_to_range_spaces[i] is not None) and (not self._dual_to_range_spaces[i].is_compatible(op.dual_to_range)):
                raise ValueError("Dual to range space not compatible")
        self._operators[key] = op

        if self._range_spaces[i] is None:
            self._range_spaces[i] = op.range

        if self._domain_spaces[j] is None:
            self._domain_spaces[j] = op.domain

        if self._dual_to_range_spaces[i] is None:
            self._dual_to_range_spaces[i] = op.dual_to_range

    def weak_form(self):

        from bempp.assembly.discrete_boundary_operator import BlockedDiscreteBoundaryOperator
        
        self.fill_complete()

        row_dimensions = np.zeros(self.ndims[0],dtype=np.int)
        column_dimensions = np.zeros(self.ndims[1],dtype=np.int)

        for i in range(self.ndims[0]):
            row_dimensions[i] = self._operators[i,0].dual_to_range.global_dof_count
        for j in range(self.ndims[1]):
            column_dimensions[j] = self._operators[0,j].domain.global_dof_count

        discrete_operator = BlockedDiscreteBoundaryOperator(
                row_dimensions,column_dimensions)

        for i in range(self.ndims[0]):
            for j in range(self.ndims[1]):
                discrete_operator[i,j] = self._operators[i,j].weak_form()
        return discrete_operator

    def __add__(self,BlockedBoundaryOperator other):

        if not self.ndims==other.ndims:
            raise ValueError("Dimensions do not match")
        
        self.fill_complete()
        other.fill_complete()

        result = BlockedBoundaryOperator(self.ndims[0],self.ndims[1])

        for i in range(self.ndims[0]):
            for j in range(self.ndims[1]):
                result[i,j] = self[i,j]+other[i,j]

        return result

    def __mul__(self, other):

        if not isinstance(self,BlockedBoundaryOperator):
            return other*self

        self.fill_complete()

        if np.isscalar(other):
            result = BlockedBoundaryOperator(self.ndims[0],self.ndims[1])
            for i in range(self.ndims[0]):
                for j in range(self.ndims[1]):
                    result[i,j] = other*self[i,j]
            return result

        if not len(other)==self.ndims[1]:
            raise NotImplementedError("Wrong type or length of multiplicator")

        result = []
        for i in range(self.ndims[0]):
            g = GridFunction(space=self.range_spaces[i],
                    coefficients=np.zeros(self.range_spaces[i].global_dof_count))
            for j in range(self.ndims[1]):
                g = g+self[i,j]*other[j]
            result.append(g)
        return result

    def __neg__(self):

        return self.__mul__(-1)

    def __sub__(self,other):

        return self.__add__(-other)







        









        






        







