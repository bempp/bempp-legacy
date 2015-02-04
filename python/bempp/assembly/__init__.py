__all__ = ['BoundaryOperatorBase','GridFunction','BlockedBoundaryOperator']
from .boundary_operator import BoundaryOperatorBase,BlockedBoundaryOperator
from .grid_function import GridFunction

__doc__="""
This package defines the basic BEM++ objects to assemble boundary operators,
functions defined on grids and domain potentiails.

Classes
-------

.. autoclass:: GridFunction
    :members: projections


"""
