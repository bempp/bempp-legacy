<%
codims = [('0','codim_zero'),('1','codim_one'),('2','codim_two')]
%>

from bempp.utils cimport Matrix
from bempp.utils.eigen cimport eigen_matrix_to_np_float64
from cython.operator cimport dereference as deref
import numpy as _np
cimport numpy as _np


% for (codim,codim_template) in codims:

cdef class Geometry${codim}:


    property corners:
        """Corners of entity"""
        def __get__(self):

            cdef Matrix[double]* c = new Matrix[double](self.dim_world,self.corner_count)
            cdef _np.ndarray res 
            self.impl_.getCorners(deref(c))
            res = eigen_matrix_to_np_float64(deref(c))
            del c
            return res

    property corner_count:
        """Number of corners of element"""
        def __get__(self):
            return self.impl_.cornerCount()

    property affine:
        """Return if element is affine"""
        def __get__(self):
            return self.impl_.affine()

    property dim:
        """" Dimension of the entity. """
        def __get__(self):
            return self.impl_.dim()

    property dim_world:
        """ Dimension of the space containing the entity. """
        def __get__(self):
            return self.impl_.dimWorld()
% endfor
