__doc__="""

This module defines the grid data structure in BEM++
and provides methods to create grids on the fly.

In order to read grids from various file formats see
the package bempp.file_interfaces.

Classes
-------

.. autoclass:: Grid

Functions
---------

.. autofunction:: grid_from_element_data
.. autofunction:: structured_grid
.. autofunction:: grid_from_sphere

"""

__all__ = ['Grid', 'structured_grid',
            'grid_from_element_data',
            'grid_from_sphere']
from .grid import Grid, structured_grid, grid_from_element_data, grid_from_sphere


