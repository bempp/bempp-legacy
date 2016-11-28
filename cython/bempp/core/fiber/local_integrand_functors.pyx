from cython.operator cimport dereference as deref

cdef extern from "bempp/core/fiber/local_integrand_functors.hpp" namespace "Fiber":
    cdef c_LocalIntegrandFunctorContainer* c_simpleTestTrialIntegrandFunctor "Fiber::simpleTestTrialIntegrandFunctor"()

cdef extern from "bempp/core/fiber/local_integrand_functors.hpp" namespace "Fiber":
    cdef c_LocalIntegrandFunctorContainer* c_maxwell3dTestTrialIntegrandFunctor "Fiber::maxwell3dTestTrialIntegrandFunctor"()

cdef extern from "bempp/core/fiber/local_integrand_functors.hpp" namespace "Fiber":
    cdef c_LocalIntegrandFunctorContainer* c_singleComponentTestTrialIntegrandFunctor "Fiber::singleComponentTestTrialIntegrandFunctor"(size_t, size_t)

cdef class LocalIntegrandFunctorContainer:

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()

def simple_test_trial_integrand_functor_ext():

    cdef LocalIntegrandFunctorContainer container = LocalIntegrandFunctorContainer()
    container.impl_.reset(c_simpleTestTrialIntegrandFunctor())
    return container


def maxwell_test_trial_integrand_functor_ext():

    cdef LocalIntegrandFunctorContainer container = LocalIntegrandFunctorContainer()
    container.impl_.reset(c_maxwell3dTestTrialIntegrandFunctor())
    return container

def single_component_test_trial_integrand_functor_ext(size_t test_component, size_t trial_component):

    cdef LocalIntegrandFunctorContainer container = LocalIntegrandFunctorContainer()
    container.impl_.reset(c_singleComponentTestTrialIntegrandFunctor(test_component, trial_component))
    return container


