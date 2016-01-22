from bempp.core.utils cimport Matrix
from bempp.core.fiber.arrays cimport _3dArray
from bempp.core.fiber.arrays cimport _4dArray

cdef extern from "bempp/fiber/basis_data.hpp" namespace "Fiber":
    cdef cppclass c_BasisData "Fiber::BasisData<double>":
        c_BasisData()
        _3dArray values
        _4dArray derivatives
        int componentCount()
        int functionCount()
        int pointCount()

cdef extern from "bempp/fiber/shapeset.hpp" namespace "Fiber":
    cdef cppclass c_Shapeset "Fiber::Shapeset<double>":
        int size()
        int order()
        void evaluate(size_t, const Matrix[double]&, int, c_BasisData&)

cdef class Shapeset:
    cdef const c_Shapeset* impl_

