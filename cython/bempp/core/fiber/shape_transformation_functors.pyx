from cython.operator cimport dereference as deref

cdef extern from "bempp/core/fiber/shape_transformation_functors.hpp" namespace "Fiber":
    cdef c_ShapeTransformationFunctorContainer* c_scalarFunctionValueFunctor "Fiber::scalarFunctionValueFunctor"()

cdef extern from "bempp/core/fiber/shape_transformation_functors.hpp" namespace "Fiber":
    cdef c_ShapeTransformationFunctorContainer* c_surfaceGrad3dFunctor "Fiber::surfaceGrad3dFunctor"()

cdef extern from "bempp/core/fiber/shape_transformation_functors.hpp" namespace "Fiber":
    cdef c_ShapeTransformationFunctorContainer* c_surfaceDiv3dFunctor "Fiber::surfaceDiv3dFunctor"()

cdef extern from "bempp/core/fiber/shape_transformation_functors.hpp" namespace "Fiber":
    cdef c_ShapeTransformationFunctorContainer* c_surfaceCurl3dFunctor "Fiber::surfaceCurl3dFunctor"()

cdef extern from "bempp/core/fiber/shape_transformation_functors.hpp" namespace "Fiber":
    cdef c_ShapeTransformationFunctorContainer* c_hdivFunctionValueFunctor "Fiber::hdivFunctionValueFunctor"()

cdef extern from "bempp/core/fiber/shape_transformation_functors.hpp" namespace "Fiber":
    cdef c_ShapeTransformationFunctorContainer* c_hcurlFunctionValueFunctor "Fiber::hcurlFunctionValueFunctor"()

cdef extern from "bempp/core/fiber/shape_transformation_functors.hpp" namespace "Fiber":
    cdef c_ShapeTransformationFunctorContainer* c_hcurlSurfaceCurlFunctor "Fiber::hcurlSurfaceCurlFunctor"()

cdef extern from "bempp/core/fiber/shape_transformation_functors.hpp" namespace "Fiber":
    cdef c_ShapeTransformationFunctorContainer* c_scalarFunctionValueTimesNormalFunctor "Fiber::scalarFunctionValueTimesNormalFunctor"()

cdef class ShapeTransformationFunctorContainer:

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()

    def argument_dimension(self):

        return deref(self.impl_).argumentDimension()

    def result_dimension(self):

        return deref(self.impl_).resultDimension(0)

def scalar_function_value_functor_ext():

    cdef ShapeTransformationFunctorContainer container = ShapeTransformationFunctorContainer()
    container.impl_.reset(c_scalarFunctionValueFunctor())
    return container

def surface_gradient_functor_ext():

    cdef ShapeTransformationFunctorContainer container = ShapeTransformationFunctorContainer()
    container.impl_.reset(c_surfaceGrad3dFunctor())
    return container

def surface_divergence_functor_ext():

    cdef ShapeTransformationFunctorContainer container = ShapeTransformationFunctorContainer()
    container.impl_.reset(c_surfaceDiv3dFunctor())
    return container

def vector_surface_curl_functor_ext():

    cdef ShapeTransformationFunctorContainer container = ShapeTransformationFunctorContainer()
    container.impl_.reset(c_surfaceCurl3dFunctor())
    return container

def hdiv_function_value_functor_ext():

    cdef ShapeTransformationFunctorContainer container = ShapeTransformationFunctorContainer()
    container.impl_.reset(c_hdivFunctionValueFunctor())
    return container

def hcurl_function_value_functor_ext():

    cdef ShapeTransformationFunctorContainer container = ShapeTransformationFunctorContainer()
    container.impl_.reset(c_hcurlFunctionValueFunctor())
    return container

def scalar_surface_curl_functor_ext():

    cdef ShapeTransformationFunctorContainer container = ShapeTransformationFunctorContainer()
    container.impl_.reset(c_hcurlSurfaceCurlFunctor())
    return container

def scalar_function_value_times_normal_functor_ext():

    cdef ShapeTransformationFunctorContainer container = ShapeTransformationFunctorContainer()
    container.impl_.reset(c_scalarFunctionValueTimesNormalFunctor())
    return container
