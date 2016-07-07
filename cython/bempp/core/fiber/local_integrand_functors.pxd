from bempp.core.utils cimport unique_ptr

cdef extern from "bempp/core/fiber/local_integrand_functors.hpp" namespace "Fiber":
    cdef cppclass c_LocalIntegrandFunctorContainer "Fiber::LocalIntegrandFunctorContainer":
        pass

cdef class LocalIntegrandFunctorContainer:
    cdef unique_ptr[c_LocalIntegrandFunctorContainer] impl_
