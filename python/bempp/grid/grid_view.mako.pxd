<%
codims = [('0','codim_zero'),('1','codim_one'),('2','codim_two')]
%>

from bempp.utils cimport unique_ptr
from bempp.grid.codim_template cimport codim_zero,codim_one,codim_two
from bempp.grid.grid cimport Grid

from bempp.grid.entity_iterator cimport c_EntityIterator
from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from bempp.utils.armadillo cimport Mat

% for (codim, codim_template) in codims:
from bempp.grid.entity_iterator cimport EntityIterator${codim}
% endfor

import numpy as np
cimport numpy as np

cdef extern from "bempp/grid/grid_view.hpp" namespace "Bempp":
    cdef cppclass c_GridView "Bempp::GridView":
        int dim() const
        int dimWorld() const
        size_t entityCount(int codim) const

        void getRawElementData(Mat[double]& vertices,
                               Mat[int]& elementCorners,
                               Mat[char]& auxData,
                               vector[int]& domainIndices)
                                    

% for (codim,codim_template) in codims:
        unique_ptr[c_EntityIterator[${codim_template}]]\
            entityIterator${codim} "entityIterator<${codim}>"() const
% endfor
        

cdef class GridView:
    cdef cbool _raw_data_is_computed 
    cdef np.ndarray _vertices
    cdef np.ndarray _elements
    cdef vector[int] _domain_indices
    cdef unique_ptr[c_GridView] impl_ 
    cdef Grid _grid
    cpdef size_t entity_count(self,int codim)
    cdef void _compute_raw_element_data(self)
 
% for (codim,codim_template) in codims:
    cpdef EntityIterator${codim} _entity_iterator${codim}(self)
% endfor
cdef GridView _grid_view_from_unique_ptr(unique_ptr[c_GridView]& c_view)
