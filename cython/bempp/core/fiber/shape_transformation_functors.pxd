from bempp.core.utils cimport unique_ptr

cdef extern from "bempp/core/fiber/shape_transformation_functors.hpp" namespace "Fiber":
    cdef cppclass c_ShapeTransformationFunctorContainer "Fiber::ShapeTransformationFunctorContainer":
        int argumentDimension()
        int resultDimension(int)

cdef class ShapeTransformationFunctorContainerExt:
    cdef unique_ptr[c_ShapeTransformationFunctorContainer] impl_
