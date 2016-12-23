"""This module contains classes and functions for shape transformations."""

class ShapeTransformationFunctorContainer(object):
    """Functor that acts on a given reference element shapeset."""

    def __init__(self, impl):
        """Constructor. Should not be called by the user."""
        self._impl = impl

    @property
    def argument_dimension(self):
        """Return the argument dimension."""
        return self._impl.argument_dimension()

    @property
    def result_dimension(self):
        """Return the result dimension."""
        return self._impl.result_dimension()

#pylint: disable=too-few-public-methods
class LocalIntegrandFunctorContainer(object):
    """Wrapper class for integrand functors of local operators."""

    def __init__(self, impl):
        """Construct functor. Should not be called by the user."""
        self._impl = impl

def scalar_function_value_functor():
    """
    Scalar function value functor.

    Identity map from scalar function values on the reference element
    to scalar function values on a grid element.

    """
    #pylint: disable=no-name-in-module
    from bempp.core.fiber.shape_transformation_functors import \
        scalar_function_value_functor_ext

    return ShapeTransformationFunctorContainer(
        scalar_function_value_functor_ext())

#pylint: disable=invalid-name
def scalar_function_value_times_normal_functor():
    """Product of scalar function value and normal vector."""
    #pylint: disable=no-name-in-module
    from bempp.core.fiber.shape_transformation_functors import \
        scalar_function_value_times_normal_functor_ext

    return ShapeTransformationFunctorContainer(
        scalar_function_value_times_normal_functor_ext())

def surface_gradient_functor():
    """Compute the surface gradient of a shapeset."""
    #pylint: disable=no-name-in-module
    from bempp.core.fiber.shape_transformation_functors import \
        surface_gradient_functor_ext

    return ShapeTransformationFunctorContainer(
        surface_gradient_functor_ext())

def surface_divergence_functor():
    """Compute the surface divergence of a shapeset."""
    #pylint: disable=no-name-in-module
    from bempp.core.fiber.shape_transformation_functors \
        import surface_divergence_functor_ext

    return ShapeTransformationFunctorContainer(
        surface_divergence_functor_ext())

def vector_surface_curl_functor():
    """Vector surface curl of a scalar function."""
    #pylint: disable=no-name-in-module
    from bempp.core.fiber.shape_transformation_functors \
        import vector_surface_curl_functor_ext

    return ShapeTransformationFunctorContainer(
        vector_surface_curl_functor_ext())

def hdiv_function_value_functor():
    """Div-conforming transformations from reference element to element."""
    #pylint: disable=no-name-in-module
    from bempp.core.fiber.shape_transformation_functors \
        import hdiv_function_value_functor_ext

    return ShapeTransformationFunctorContainer(
        hdiv_function_value_functor_ext())

def hcurl_function_value_functor():
    """Curl conforming transformations from reference element to element."""
    #pylint: disable=no-name-in-module
    from bempp.core.fiber.shape_transformation_functors \
        import hcurl_function_value_functor_ext

    return ShapeTransformationFunctorContainer(
        hcurl_function_value_functor_ext())

def scalar_surface_curl_functor():
    """Scalar surface curl of a vector function."""
    #pylint: disable=no-name-in-module
    from bempp.core.fiber.shape_transformation_functors import \
        scalar_surface_curl_functor_ext

    return ShapeTransformationFunctorContainer(
        scalar_surface_curl_functor_ext())

def simple_test_trial_integrand_functor():
    """Simple integral between test and trial functions."""
    #pylint: disable=no-name-in-module
    from bempp.core.fiber.local_integrand_functors import \
        simple_test_trial_integrand_functor_ext

    return LocalIntegrandFunctorContainer(
        simple_test_trial_integrand_functor_ext())

def single_component_test_trial_integrand_functor(
        test_component, trial_component):
    """
    Multiply one component of the trial with one component of test function.

    Parameters
    ----------
    test_component : int
        index of the component of the test function.
    trial_component : int
        index of the component of the trial function.

    """
    #pylint: disable=no-name-in-module
    from bempp.core.fiber.local_integrand_functors import \
        single_component_test_trial_integrand_functor_ext

    return LocalIntegrandFunctorContainer(
        single_component_test_trial_integrand_functor_ext(
            test_component, trial_component))

def maxwell_test_trial_integrand_functor():
    """
    Multiplication of test and trial component in twisted inner prodct.

    If f is the trial function and g the test function this functor
    evaluates the inner product < f x n, g>, where n is the normal
    vector on the element. This is the standard inner product for
    Maxwell computations.

    """
    #pylint: disable=no-name-in-module
    from bempp.core.fiber.local_integrand_functors import \
        maxwell_test_trial_integrand_functor_ext

    return LocalIntegrandFunctorContainer(
        maxwell_test_trial_integrand_functor_ext())
