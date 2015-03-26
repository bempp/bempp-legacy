<%
codims = [('0','codim_zero'),('1','codim_one'),('2','codim_two')]
%>

from libcpp cimport bool as cbool
from bempp.utils.eigen cimport Matrix

cdef extern from "bempp/grid/geometry.hpp" namespace "Bempp":
    cdef cppclass c_Geometry "Bempp::Geometry":
        int dim() const
        int dimWorld() const
        cbool affine() const
        int cornerCount() const
        void getCorners(Matrix[double]& c) const

% for (codim,codim_template) in codims:

from bempp.grid.entity cimport Entity${codim}

cdef class Geometry${codim}:
    cdef const c_Geometry * impl_
    cdef Entity${codim} _entity

% endfor
