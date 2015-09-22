from bempp.core.utils cimport unique_ptr
from bempp.core.grid.codim_template cimport codim_zero,codim_one,codim_two
from bempp.core.grid.grid cimport Grid
from bempp.core.grid.index_set cimport c_IndexSet, IndexSet

from bempp.core.grid.entity_iterator cimport c_EntityIterator
from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from bempp.core.utils cimport Matrix, Vector

from bempp.core.grid.entity_iterator cimport EntityIterator0
from bempp.core.grid.entity_iterator cimport EntityIterator1
from bempp.core.grid.entity_iterator cimport EntityIterator2

import numpy as np
cimport numpy as np

cdef extern from "bempp/grid/grid_view.hpp" namespace "Bempp":
    cdef cppclass c_GridView "Bempp::GridView":
        int dim() const
        int dimWorld() const
        size_t entityCount(int codim) const
        const c_IndexSet& indexSet() const

        void getRawElementData(Matrix[double]& vertices,
                               Matrix[int]& elementCorners,
                               Matrix[char]& auxData,
                               vector[int]& domainIndices)

        unique_ptr[c_EntityIterator[codim_zero]]\
            entityIterator0 "entityIterator<0>"() const
        unique_ptr[c_EntityIterator[codim_one]]\
            entityIterator1 "entityIterator<1>"() const
        unique_ptr[c_EntityIterator[codim_two]]\
            entityIterator2 "entityIterator<2>"() const
        

cdef class GridView:
    cdef cbool _raw_data_is_computed 
    cdef np.ndarray _vertices
    cdef np.ndarray _elements
    cdef vector[int] _domain_indices
    cdef unique_ptr[c_GridView] impl_ 
    cdef Grid _grid
    cpdef size_t entity_count(self,int codim)
    cpdef IndexSet index_set(self)
    cdef void _compute_raw_element_data(self)
    cpdef EntityIterator0 _entity_iterator0(self)
    cpdef EntityIterator1 _entity_iterator1(self)
    cpdef EntityIterator2 _entity_iterator2(self)

cdef GridView _grid_view_from_unique_ptr(unique_ptr[c_GridView]& c_view)
