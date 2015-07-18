__all__ = ['GridFunction',
            'BlockedBoundaryOperator',
            'InverseSparseDiscreteBoundaryOperator']
from .boundary_operator import BlockedBoundaryOperator
from .grid_function import GridFunction
from .discrete_boundary_operator import InverseSparseDiscreteBoundaryOperator

__doc__="""
This package defines the basic BEM++ objects to assemble boundary operators,
functions defined on grids and domain potentiails.

Classes
-------

.. autoclass:: GridFunction
    :members: projections


"""
