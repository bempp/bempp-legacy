<%
codims = [('0','codim_zero'),('1','codim_one'),('2','codim_two')]
%>

from bempp.utils.armadillo cimport Mat
from cython.operator cimport dereference as deref
import numpy as np
cimport numpy as np


% for (codim,codim_template) in codims:

cdef class Geometry${codim}:


    property corners:
        """Corners of entity"""
        def __get__(self):

            cdef np.ndarray c_np = np.empty((self.dim_world,self.corner_count),dtype='float64',order='F')
            cdef double[::1,:] c_np_view = c_np
            cdef Mat[double]* c = new Mat[double](&c_np_view[0,0],self.dim_world,self.corner_count,False,True)
            self.impl_.getCorners(deref(c))
            del c
            return c_np

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
