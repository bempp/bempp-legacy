""" Boundary Element Method package BEM++ """
__all__ = [
    'Grid', 'config', '__version__', 'Options', 'space', 'scalar_spaces'
]

from grid import Grid
from . import space
from .space import scalar as scalar_spaces
import config
from .options import Options
# A fair number of packages expose version info via similar variable
from config import version as __version__
