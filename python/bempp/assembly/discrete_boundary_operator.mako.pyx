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
% for pyvalue in dtypes:
from bempp.utils.armadillo cimport armadillo_to_np_${pyvalue}
% endfor
cimport numpy as np
import numpy as np
cimport cython


cdef class DiscreteBoundaryOperatorBase:

    property dtype:
        def __get__(self):
            return self._value_type

% for pyvalue,cyvalue in dtypes.items():
    cdef void _apply_${pyvalue}(self,
            TranspositionMode trans, 
            np.ndarray[${scalar_cython_type(cyvalue)},ndim=2,mode='fortran'] x_in, 
            np.ndarray[${scalar_cython_type(cyvalue)},ndim=2,mode='fortran'] y_inout, 
            ${scalar_cython_type(cyvalue)} alpha,
            ${scalar_cython_type(cyvalue)} beta):

        raise NotImplementedError("Method _apply_${pyvalue} is not implemented.")
% endfor

    cpdef np.ndarray as_matrix(self):

        return self._as_matrix()

    cpdef np.ndarray _as_matrix(self):

        raise NotImplementedError("Method not implemented.")


    cpdef object apply(self,np.ndarray x,np.ndarray y,object transpose,object alpha, object beta):

        if not (x.dtype==self.dtype and y.dtype==self.dtype):
            raise ValueError("Wrong dtype of input arrays")

        if not x.flags['F_CONTIGUOUS'] or not y.flags['F_CONTIGUOUS']:
            raise ValueError("Input arrays must be in Fortran order")

% for pyvalue,cyvalue in dtypes.items():

        if self._value_type == "${pyvalue}":
            self._apply_${pyvalue}(transposition_mode(convert_to_bytes(transpose)),
                    x,y,alpha,beta)
            return None
% endfor


    def matvec(self,np.ndarray x):

        cdef np.ndarray x_in
        cdef np.ndarray y_inout

        cdef int rows = self.shape[0]
        cdef int cols = self.shape[1]

        if (x.ndim==1):
            x_in = x.reshape((-1,1),order='F').astype(self.dtype,
                    order='F',casting='safe',copy=False)
        elif (x.ndim==2):
            x_in = x.astype(self.dtype,order='F',casting='safe',copy=False)
        else:
            raise ValueError('x must have at most two dimensions')

        y = np.zeros((rows,x.shape[1]),dtype=self.dtype,order='F')

        self.apply(x_in,y,'no_transpose',1.0,0.0)
        return y



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
            if self._value_type=="${pyvalue}":
                rows = deref(self._impl_${pyvalue}_).rowCount()
                cols = deref(self._impl_${pyvalue}_).columnCount()
                return (rows,cols)
% endfor
            raise ValueError("Unknown value type")
    
    cpdef np.ndarray _as_matrix(self):

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
    cdef DiscreteBoundaryOperatorBase op

% for pyvalue,cyvalue in dtypes.items():
    cdef ${scalar_cython_type(cyvalue)} alpha_${pyvalue}
% endfor    

    def __cinit__(self,DiscreteBoundaryOperatorBase op,object alpha):

        self.op = op
        self._value_type = op._value_type
% for pyvalue in dtypes:
        if self._value_type=="${pyvalue}":
            self.alpha_${pyvalue} = alpha
% endfor


    def __init__(self,object op, object alpha):

        pass

    cpdef np.ndarray _as_matrix(self):

% for pyvalue in dtypes:
        if self._value_type=="${pyvalue}":
            return self.alpha_${pyvalue}*self.op.as_matrix()
% endfor 
        raise ValueError("Object has unknown dtype")

    property shape:

        def __get__(self):
            return self.op.shape

% for pyvalue,cyvalue in dtypes.items():
    cdef void _apply_${pyvalue}(self,
            TranspositionMode trans, 
            np.ndarray[${scalar_cython_type(cyvalue)},ndim=2,mode='fortran'] x_in, 
            np.ndarray[${scalar_cython_type(cyvalue)},ndim=2,mode='fortran'] y_inout, 
            ${scalar_cython_type(cyvalue)} alpha,
            ${scalar_cython_type(cyvalue)} beta):


        self.op._apply_${pyvalue}(trans,x_in,y_inout,self.alpha_${pyvalue}*alpha,beta)
% endfor





