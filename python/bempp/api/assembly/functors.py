
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


def scalarFunctionValueFunctor():

    from bempp.core.fiber.shape_transformation_functors import ShapeTransformationFunctorContainerExt
    from bempp.core.fiber.shape_transformation_functors import scalarFunctionValueFunctorExt

    return ShapeTransformationFunctorContainer(
            scalarFunctionValueFunctorExt())



def simpleTestTrialIntegrandFunctor():

    from bempp.core.fiber.local_integrand_functors import LocalIntegrandFunctorContainerExt
    from bempp.core.fiber.local_integrand_functors import simpleTestTrialIntegrandFunctorExt

    return LocalIntegrandFunctorContainer(
            simpleTestTrialIntegrandFunctorExt)



