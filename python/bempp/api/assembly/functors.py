
"""This module contains classes and functions for shape transformation functors."""

class ShapeTransformationFunctorContainer(object):
    """Wrapper class for a functor that acts on a given reference element shapeset."""


    def __init__(self, impl):

        self._impl = impl

    @property
    def argumentDimension(self):
        """Return the argument dimension."""

        return self._impl.argumentDimension()

    @property
    def resultDimension(self):
        """Return the result dimension."""

        return self._impl.resultDimension()

class LocalIntegrandFunctorContainer(object):
    """Wrapper class for integrand functors of local operators."""

    def __init__(self, impl):

        self._impl = impl

def scalar_function_value_functor():

    from bempp.core.fiber.shape_transformation_functors import scalar_function_value_functor_ext

    return ShapeTransformationFunctorContainer(
            scalar_function_value_functor_ext())

def scalar_function_value_times_normal_functor():

    from bempp.core.fiber.shape_transformation_functors import scalar_function_value_times_normal_functor_ext

    return ShapeTransformationFunctorContainer(
            scalar_function_value_times_normal_functor_ext())

def surface_gradient_functor():

    from bempp.core.fiber.shape_transformation_functors import surface_gradient_functor_ext

    return ShapeTransformationFunctorContainer(
            surface_gradient_functor_ext())

def surface_divergence_functor():

    from bempp.core.fiber.shape_transformation_functors import surface_divergence_functor_ext

    return ShapeTransformationFunctorContainer(
            surface_divergence_functor_ext())

def vector_surface_curl_functor():

    from bempp.core.fiber.shape_transformation_functors import vector_surface_curl_functor_ext

    return ShapeTransformationFunctorContainer(
            vector_surface_curl_functor_ext())

def hdiv_function_value_functor():

    from bempp.core.fiber.shape_transformation_functors import hdiv_function_value_functor_ext

    return ShapeTransformationFunctorContainer(
            hdiv_function_value_functor_ext())

def hcurl_function_value_functor():

    from bempp.core.fiber.shape_transformation_functors import hcurl_function_value_functor_ext

    return ShapeTransformationFunctorContainer(
            hcurl_function_value_functor_ext())

def scalar_surface_curl_functor():

    from bempp.core.fiber.shape_transformation_functors import scalar_surface_curl_functor_ext

    return ShapeTransformationFunctorContainer(
            scalar_surface_curl_functor_ext())

def simple_test_trial_integrand_functor():

    from bempp.core.fiber.local_integrand_functors import simple_test_trial_integrand_functor_ext

    return LocalIntegrandFunctorContainer(
            simple_test_trial_integrand_functor_ext())

def single_component_test_trial_integrand_functor(test_component, trial_component):

    from bempp.core.fiber.local_integrand_functors import single_component_test_trial_integrand_functor_ext

    return LocalIntegrandFunctorContainer(
            single_component_test_trial_integrand_functor_ext(test_component, trial_component))

def maxwell_test_trial_integrand_functor():

    from bempp.core.fiber.local_integrand_functors import maxwell_test_trial_integrand_functor_ext

    return LocalIntegrandFunctorContainer(
            maxwell_test_trial_integrand_functor_ext())


