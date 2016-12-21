"""The grid package implements interfaces to the underlying grid structures."""

from bempp.api.grid.grid import structured_grid, grid_from_element_data
from bempp.api.grid.grid_factory import GridFactory


__all__ = ['structured_grid',
           'grid_from_element_data',
           'GridFactory']
