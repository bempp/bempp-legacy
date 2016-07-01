from cython.operator cimport dereference as deref

cdef extern from "bempp/core/fiber/local_integrand_functors.hpp" namespace "Fiber":
    cdef c_LocalIntegrandFunctorContainer* c_simpleTestTrialIntegrandFunctor "Fiber::simpleTestTrialIntegrandFunctor"()

cdef class LocalIntegrandFunctorContainerExt:

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()

def simpleTestTrialIntegrandFunctorExt():

    cdef LocalIntegrandFunctorContainerExt container = LocalIntegrandFunctorContainerExt()
    container.impl_.reset(c_simpleTestTrialIntegrandFunctor())
    return container
