from cython.operator cimport dereference as deref

cdef extern from "bempp/core/fiber/shape_transformation_functors.hpp" namespace "Fiber":
    cdef c_ShapeTransformationFunctorContainer* c_scalarFunctionValueFunctor "Fiber::scalarFunctionValueFunctor"()

cdef class ShapeTransformationFunctorContainerExt:

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()

    def argumentDimension(self):

        return deref(self.impl_).argumentDimension()

    def resultDimension(self):

        return deref(self.impl_).resultDimension(0)

def scalarFunctionValueFunctorExt():

    cdef ShapeTransformationFunctorContainerExt container = ShapeTransformationFunctorContainerExt()
    container.impl_.reset(c_scalarFunctionValueFunctor())
    return container




