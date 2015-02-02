<%
from data_types import dtypes, scalar_cython_type
%>

from bempp.utils.armadillo cimport Mat
from bempp.utils.enum_types cimport transposition_mode
from bempp.utils cimport complex_float,complex_double
from cython.operator cimport dereference as deref
from bempp.utils.enum_types cimport TranspositionMode
cimport bempp.utils.enum_types as enums
from bempp.utils.byte_conversion import convert_to_bytes
from bempp.utils cimport shared_ptr, static_pointer_cast
% for pyvalue in dtypes:
from bempp.utils.armadillo cimport armadillo_to_np_${pyvalue}
% endfor
from bempp.utils import combined_type
cimport numpy as np
import numpy as np
cimport cython
from libcpp cimport bool


cdef class DiscreteBoundaryOperatorBase:

    property dtype:
        def __get__(self):
            return self._dtype

    property shape:
        def __get__(self):
            raise NotImplementedError("Method not implemented.")

    def as_matrix(self):

        raise NotImplementedError("Method not implemented.")

    def matvec(self,np.ndarray x):

        raise NotImplementedError("Method not implemented.")

    def matmat(self,np.ndarray x):

        return self.matvec(x)

    def __call__(self,np.ndarray x):

        return self.matvec(x)

    def __mul__(self,object x):

        if not isinstance(self,DiscreteBoundaryOperatorBase):
            return x*self

        if isinstance(x,DiscreteBoundaryOperatorBase):
            return _ProductDiscreteBoundaryOperator(self,x)
        elif np.isscalar(x):
            return _ScaledDiscreteBoundaryOperator(self,x)
        else:
            return self.matvec(x)

    def dot(self,object other):
        
        return self.__mul__(other)

    def __add__(self,DiscreteBoundaryOperatorBase other):

        return _SumDiscreteBoundaryOperator(self,other)

    def __neg__(self):

        return _ScaledDiscreteBoundaryOperator(self,-1.0)

    def __sub__(self,DiscreteBoundaryOperatorBase x):
        return self.__add__(-x)
    
    def __repr__(self):
        
        M,N = self.shape
        dt = 'dtype=' + str(self.dtype)
        return '<%text><%dx%d %s with %s></%text>' % (M, N, self.__class__.__name__, dt)

cdef class ZeroDiscreteBoundaryOperator(DiscreteBoundaryOperatorBase):

    cdef object _shape

    def __cinit__(self,int M, int N):
        pass

    def __init__(self,int M, int N):
        self._shape = (M,N)
        self._dtype = np.dtype('float64')

    property shape:
        def __get__(self):
            return self._shape

    def as_matrix(self):

        return np.zeros(self.shape,dtype=self.dtype)

    def matvec(self, np.ndarray x):

        cdef bool is_reshaped=False
        if not x.shape[0]==self.shape[1]:
            return ValueError("Wrong dimensions.")

        if x.ndim==1:
            x = x.reshape((-1,1))
            is_reshaped=True

        cdef np.ndarray result = np.zeros((self.shape[0],x.shape[1]),
                dtype=x.dtype)
        if is_reshaped:
            result = result.ravel()
        return result



cdef class DiscreteBoundaryOperator(DiscreteBoundaryOperatorBase):

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    property shape:

        def __get__(self):
            cdef unsigned int rows 
            cdef unsigned int cols

% for pyvalue,cyvalue in dtypes.items():
            if self.dtype=="${pyvalue}":
                rows = deref(self._impl_${pyvalue}_).rowCount()
                cols = deref(self._impl_${pyvalue}_).columnCount()
                return (rows,cols)
% endfor
            raise ValueError("Unknown value type")
    
    def as_matrix(self):

% for pyvalue in dtypes:

        if self.dtype=="${pyvalue}":
            return self._as_matrix_${pyvalue}()
% endfor

        raise ValueError("Unknown value type")


% for pyvalue,cyvalue in dtypes.items():
    cdef np.ndarray _as_matrix_${pyvalue}(self):

        cdef Mat[${cyvalue}] mat_data = deref(self._impl_${pyvalue}_).asMatrix()
        return armadillo_to_np_${pyvalue}(mat_data)

% endfor

    def _apply(self,np.ndarray x,np.ndarray y,object transpose,object alpha, object beta):

        if not (x.dtype==self.dtype and y.dtype==self.dtype):
            raise ValueError("Wrong dtype of input arrays")

        if not x.flags['F_CONTIGUOUS'] or not y.flags['F_CONTIGUOUS']:
            raise ValueError("Input arrays must be in Fortran order")

% for pyvalue,cyvalue in dtypes.items():

        if self.dtype == "${pyvalue}":
            self._apply_${pyvalue}(transposition_mode(convert_to_bytes(transpose)),
                    x,y,alpha,beta)
            return None
% endfor


    def matvec(self,np.ndarray x):

        if self.dtype=='float64' and np.iscomplexobj(x):
            return self*np.real(x)+1j*(self*np.imag(x))

        cdef np.ndarray x_in
        cdef np.ndarray y_inout

        cdef int rows = self.shape[0]
        cdef int cols = self.shape[1]

        cdef bool is_reshaped = False

        if (x.ndim==1):
            x_in = x.reshape((-1,1),order='F').astype(self.dtype,
                    order='F',casting='safe',copy=False)
            is_reshaped = True
        elif (x.ndim==2):
            x_in = x.astype(self.dtype,order='F',casting='safe',copy=False)
        else:
            raise ValueError('x must have at most two dimensions')

        y = np.zeros((rows,x_in.shape[1]),dtype=self.dtype,order='F')

        self._apply(x_in,y,'no_transpose',1.0,0.0)
        if is_reshaped:
            y = y.ravel()
        return y
        
% for pyvalue,cyvalue in dtypes.items():
    cdef void _apply_${pyvalue}(self,
            TranspositionMode trans, 
            np.ndarray[${scalar_cython_type(cyvalue)},ndim=2,mode='fortran'] x_in, 
            np.ndarray[${scalar_cython_type(cyvalue)},ndim=2,mode='fortran'] y_inout, 
            ${scalar_cython_type(cyvalue)} alpha,
            ${scalar_cython_type(cyvalue)} beta):


        cdef int rows = self.shape[0]
        cdef int cols = self.shape[1]
        cdef int xrows = x_in.shape[0]
        cdef int xcols = x_in.shape[1]
        cdef int yrows = y_inout.shape[0]
        cdef int ycols = y_inout.shape[1]

        cdef Mat[${cyvalue}]* arma_${pyvalue}_buff_x
        cdef Mat[${cyvalue}]* arma_${pyvalue}_buff_y

        if (trans==enums.no_transpose or trans==enums.conjugate):

            if not (rows==yrows and cols ==xrows):
                raise ValueError("Wrong dimensions")

            if not (xcols==ycols):
                raise ValueError("Wrong dimensions")
        elif (trans==enums.transpose or trans==enums.conjugate_transpose):

            if not (cols==yrows and rows ==xrows):
                raise ValueError("Wrong dimensions")

            if not (xcols==ycols):
                raise ValueError("Wrong dimensions")
        else:
            raise ValueError("Unknown transposition mode")


% if pyvalue=='complex128' or pyvalue=='complex64':
        cdef ${cyvalue} cpp_alpha = ${cyvalue}(alpha.real,alpha.imag)
        cdef ${cyvalue} cpp_beta = ${cyvalue}(beta.real,beta.imag)
% else:
        cdef ${cyvalue} cpp_alpha = alpha
        cdef ${cyvalue} cpp_beta = beta
% endif

        arma_${pyvalue}_buff_x = new Mat[${cyvalue}](<${cyvalue}*>&x_in[0,0],xrows,xcols,False,True)
        arma_${pyvalue}_buff_y = new Mat[${cyvalue}](<${cyvalue}*>&y_inout[0,0],yrows,ycols,False,True)

        deref(self._impl_${pyvalue}_).apply(trans,
                deref(arma_${pyvalue}_buff_x),
                deref(arma_${pyvalue}_buff_y),
                cpp_alpha,cpp_beta)

        del arma_${pyvalue}_buff_y
        del arma_${pyvalue}_buff_x
% endfor


cdef class _ScaledDiscreteBoundaryOperator(DiscreteBoundaryOperatorBase):
    cdef DiscreteBoundaryOperatorBase _op
    cdef object _alpha

    def __cinit__(self,DiscreteBoundaryOperatorBase op,object alpha):
        pass

    def __init__(self,DiscreteBoundaryOperatorBase op, object alpha):

        self._op = op
        self._alpha = 1.0*alpha # make sure it is not integer 
        self._dtype = combined_type(np.dtype(type(self._alpha)),op.dtype)

    def as_matrix(self):

        return self._alpha*self._op.as_matrix()

    property shape:

        def __get__(self):
            return self._op.shape

    def matvec(self,np.ndarray x):

        return self._alpha*(self._op*x)

cdef class _SumDiscreteBoundaryOperator(DiscreteBoundaryOperatorBase):

    cdef DiscreteBoundaryOperatorBase _op1
    cdef DiscreteBoundaryOperatorBase _op2

    def __cinit__(self,DiscreteBoundaryOperatorBase op1, DiscreteBoundaryOperatorBase op2):
        pass

    def __init__(self,DiscreteBoundaryOperatorBase op1, DiscreteBoundaryOperatorBase op2):
        if not op1.shape == op2.shape:
            raise ValueError("Both operators must have the same shape")

        self._op1 = op1
        self._op2 = op2

        self._dtype = combined_type(self._op1.dtype,self._op2.dtype)

    property shape:

        def __get__(self):

            return self._op1.shape

    def as_matrix(self):

        return self._op1.as_matrix()+self._op2.as_matrix()

    def matvec(self,np.ndarray x):

        return self._op1*x+self._op2*x


cdef class _ProductDiscreteBoundaryOperator(DiscreteBoundaryOperatorBase):

    cdef DiscreteBoundaryOperatorBase _op1
    cdef DiscreteBoundaryOperatorBase _op2

    def __cinit__(self,DiscreteBoundaryOperatorBase op1, DiscreteBoundaryOperatorBase op2):
        pass

    def __init__(self,DiscreteBoundaryOperatorBase op1, DiscreteBoundaryOperatorBase op2):
        if not op1.shape[1]==op2.shape[0]:
            raise ValueError("Incompatible Dimensions.")

        self._op1 = op1
        self._op2 = op2

        self._dtype = combined_type(self._op1.dtype,self._op2.dtype)

    def as_matrix(self):
        return self._op1.as_matrix()*self._op2.as_matrix()

    def matvec(self,np.ndarray x):

        return self._op1*(self._op2*x)

    property shape:

        def __get__(self):
            return (self._op1.shape[0],self._op2.shape[1])

cdef class SparseDiscreteBoundaryOperator(DiscreteBoundaryOperatorBase):

    def __cinit__(self,op):
        pass

    def __init__(self,op):
        from scipy.sparse import csr_matrix

        if not isinstance(op,csr_matrix):
            raise ValueError("op must be of type scipy.sparse.csr.csr_matrix")

        self._op = op
        self._dtype = self._op.dtype

    def as_matrix(self):

        return self._op.todense()

    def matvec(self,x):

        return self._op*x

    def __add__(self,DiscreteBoundaryOperatorBase other):

        if isinstance(other,SparseDiscreteBoundaryOperator):
            return SparseDiscreteBoundaryOperator(self.sparse_operator+other.sparse_operator)

        else:
            return super(SparseDiscreteBoundaryOperator,self).__add__(other)

    def __mul__(self,object x):

        if not isinstance(self,SparseDiscreteBoundaryOperator):
            return x*self

        if np.isscalar(x):
            return SparseDiscreteBoundaryOperator(x*self.sparse_operator)
        return super(SparseDiscreteBoundaryOperator,self).__mul__(x)

    def __neg__(self):

        return SparseDiscreteBoundaryOperator(-self.sparse_operator)

    def __sub__(self,other):

        if isinstance(other,SparseDiscreteBoundaryOperator):
            return SparseDiscreteBoundaryOperator(self.sparse_operator-other.sparse_operator)

        else:
            return super(SparseDiscreteBoundaryOperator,self).__sub__(other)


    property sparse_operator:
        """ The SciPy sparse matrix representation of the operator """

        def __get__(self):
            return self._op

    property shape:

        def __get__(self):
            return self._op.shape

    property dtype:
        def __get__(self):
            return self._op.dtype

cdef class DenseDiscreteBoundaryOperator(DiscreteBoundaryOperator):

    def __cinit__(self):
        pass

    def __init__(self):
        self._array_view = None

    cdef object _init_array_view(self):
        """ Initialize the view on the dense operator via a Numpy Array """
        if self._array_view is not None:
            return
% for pyvalue,cyvalue in dtypes.items():
        if self.dtype=="${pyvalue}":
            self._array_view = py_array_from_dense_operator[${cyvalue}](self._impl_${pyvalue}_)
            return
% endfor

        raise ValueError("Unknown data type")

    property numpy_view:

        def __get__(self):
            self._init_array_view()
            return self._array_view[:]

cdef class BlockedDiscereteBoundaryOperator(DiscreteBoundaryOperatorBase):

    cdef np.ndarray _row_dimensions
    cdef np.ndarray _column_dimensions

    cdef np.ndarray _row_sums
    cdef np.ndarray _column_sums

    cdef np.ndarray _operators

    cdef object _shape

    def __cinit__(self,np.ndarray[int,ndim=1,mode='fortran'] row_dimensions, 
            np.ndarray[int,ndim=1,mode='fortran'] column_dimensions):
        pass

    def __init__(self,np.ndarray[int,ndim=1,mode='fortran'] row_dimensions, 
            np.ndarray[int,ndim=1,mode='fortran'] column_dimensions):
        self._row_dimensions = row_dimensions
        self._column_dimensions = column_dimensions
        self._row_sums = np.hstack(([0],np.cumsum(self._row_dimensions)))
        self._column_sums = np.hstack(([0],np.cumsm(self._column_dimensions)))

        self._shape = (len(row_dimensions),len(column_dimensions))

        self._operators = np.empty(self._shape,dtype=np.object)

        for row,i in enumerate(self._row_dimensions):
            for col,j in enumerate(self._column_dimensions):
                self._operators[i,j] = ZeroDiscreteBoundaryOperator(
                        row,col)
                
    property shape:

        def __get__(self):
            return self._shape

    property dtype:

        def __get__(self):
            cdef object dt = np.dtype('float64')
            for op in self._operators.ravel():
                dt = combined_type(dt,op.dtype)
            return dt

    def __get__item(self,key):
        return self._operators[key]

    def __setitem(self,key, DiscreteBoundaryOperatorBase op):

        if not (op.shape==self._operators[key].shape):
            raise ValueError("Wrong dimensions")

        self._operators[key] = op

    def as_matrix(self):

        res = np.empty(self.shape,dtype=self.dtype)
        for row_dim,i in enumerate(self._row_dimensions):
            for col_dim,j in enumerate(self._column_dimensions):
                res[self._row_sums[i]:self._row_sums[i+1],
                        self._column_sums[i]:self._column_sums[i+1]]=self._operators[i,j].as_matrix()
        return res


    def matvec(self,np.ndarray x):

        res = np.zeros((self.shape[0],x.shape[1]),dtype=self.dtype)

        for i in range(len(self._row_dimensions)):
            for j in range(len(self._column_dimensions)):
                res[self._row_sums[i]:self._row_sums[i+1],:] = res[self._row_sums[i]:self._row_sums[i+1],:] + self[i,j]*x[self._column_sums[j]:self._column_sums[j+1],:]
        return res





    
        
