<%
from data_types import dtypes, scalar_cython_type
%>

from bempp.utils.armadillo cimport Mat
from bempp.utils.enum_types cimport transposition_mode
from bempp.utils cimport complex_float,complex_double
from cython.operator cimport dereference as deref
from bempp.utils.enum_types cimport TranspositionMode

cdef class DiscreteBoundaryOperatorBase:

% for pyvalue,cyvalue in dtypes.items():
    cdef void _apply_${pyvalue}(self,
            TranspositionMode trans, const Mat[${cyvalue}]& x_in,
            Mat[${cyvalue}]& y_inout, ${scalar_cython_type(cyvalue)} alpha,
            ${scalar_cython_type(cyvalue)} beta):

        raise NotImplementedError("Method _apply_${pyvalue} is not implemented.")
% endfor

    cpdef object apply(self,np.ndarray x,np.ndarray y,object transpose,object alpha, object beta):

        cdef int rows = self.shape[0]
        cdef int cols = self.shape[1]
        cdef int xrows = x.shape[0]
        cdef int xcols = x.shape[1]
        cdef int yrows = y.shape[0]
        cdef int ycols = y.shape[1]

% for pyvalue,cyvalue in dtypes.items():
        cdef np.ndarray[${scalar_cython_type(cyvalue)},ndim=2,mode='fortran'] ${pyvalue}_buff_x = x
        cdef np.ndarray[${scalar_cython_type(cyvalue)},ndim=2,mode='fortran'] ${pyvalue}_buff_y = y
        cdef Mat[${cyvalue}]* arma_${pyvalue}_buff_x
        cdef Mat[${cyvalue}]* arma_${pyvalue}_buff_y
% endfor        

        if not (rows==y.shape[0] and cols ==x.shape[0]):
            raise ValueError("Wrong dimensions")

        if not (x.shape[1]==y.shape[1]):
            raise ValueError("Wrong dimensions")

% for pyvalue,cyvalue in dtypes.items():

        if self._value_type == "${pyvalue}":
            arma_${pyvalue}_buff_x = new Mat[${cyvalue}](<${cyvalue}*>&${pyvalue}_buff_x[0,0],xrows,xcols,False,True)
            arma_${pyvalue}_buff_y = new Mat[${cyvalue}](<${cyvalue}*>&${pyvalue}_buff_y[0,0],yrows,ycols,False,True)
            self._apply_${pyvalue}(transposition_mode(transpose),x,y,alpha,beta)
            del arma_${pyvalue}_buff_y
            del arma_${pyvalue}_buff_x

            return None
% endfor

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

% for pyvalue,cyvalue in dtypes.items():
    cdef void _apply_${pyvalue}(self,
            TranspositionMode trans, const Mat[${cyvalue}]& x_in,
            Mat[${cyvalue}]& y_inout, 
            ${scalar_cython_type(cyvalue)} alpha,
            ${scalar_cython_type(cyvalue)} beta):

% if pyvalue=='complex128' or pyvalue=='complex64':
        cdef ${cyvalue} cpp_alpha = ${cyvalue}(alpha.real,alpha.imag)
        cdef ${cyvalue} cpp_beta = ${cyvalue}(beta.real,beta.imag)
% else:
        cdef ${cyvalue} cpp_alpha = alpha
        cdef ${cyvalue} cpp_beta = beta
% endif
        deref(self._impl_${pyvalue}_).apply(trans,x_in,y_inout,cpp_alpha,cpp_beta)

% endfor
