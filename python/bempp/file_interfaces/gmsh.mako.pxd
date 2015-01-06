from bempp.utils cimport shared_ptr
from bempp.grid.grid cimport c_Grid
from libcpp.vector cimport vector
from libcpp.string cimport string
from bempp.utils cimport catch_exception


cdef extern from "bempp/io/gmsh.hpp" namespace "Bempp":
    cdef cppclass c_GmshIo "Bempp::GmshIo":
        c_GmshIo(shared_ptr[const c_Grid]& grid) except+catch_exception
        c_GmshIo(string fileName, int physicalEntity) except+catch_exception

        shared_ptr[const c_Grid] grid() const
        vector[int] nodePermutation() const
        vector[int] elementPermutation() const
        vector[int] inverseNodePermutation() const
        vector[int] inverseElementPermutation() const

        void write(string fileName) except+catch_exception

cdef class Gmsh:
    cdef shared_ptr[c_GmshIo] impl_








        
