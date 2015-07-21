from bempp.utils cimport catch_exception
from bempp.utils cimport shared_ptr
from bempp.grid.grid cimport c_Grid
from bempp.grid.grid cimport Grid
from libcpp cimport bool as cbool

cdef extern from "bempp/grid/grid_factory.hpp" namespace "Bempp":
    cdef cppclass c_GridFactory "Bempp::GridFactory":
        c_GridFactory()
        void insertVertex(double x, double y, double z) except +catch_exception
        void insertElement(unsigned int v0, unsigned int v1, unsigned int v2, int domain_index) except +catch_exception
        shared_ptr[c_Grid] finalize() except +catch_exception


cdef class GridFactory(object):

    cdef c_GridFactory _grid_factory
    cdef cbool finalized

    def __cinit__(self):
        self.finalized = False

    def __init__(self):
        pass

    def insert_vertex(self, vertex):
        if self.finalized:
            raise Exception("Grid is already finalized.")
        self._grid_factory.insertVertex(vertex[0],vertex[1],vertex[2])

    def insert_element(self, element, domain_index=0):
        if self.finalized:
            raise Exception("Grid is already finalized.")
        self._grid_factory.insertElement(element[0],element[1],element[2],domain_index)

    def finalize(self):
        if self.finalized:
            raise Exception("Grid is already finalized.")
        cdef Grid grid = Grid.__new__(Grid)
        grid.impl_ = self._grid_factory.finalize()
        return grid



    




