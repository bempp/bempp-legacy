""" Boundary Element Method package BEM++ """
from __future__ import print_function

__all__ = ['grid_from_sphere',
           'grid_from_element_data',
           'structured_grid',
           'function_space',
           'GridFunction',
           'global_parameters',
           'BlockedBoundaryOperator',
           'ZeroBoundaryOperator',
           'export',
           'GridFactory',
           'import_grid',
           'operators',
           'shapes',
           'applications',
           'FileReader',
           'generate_block_cluster_tree',
           'InverseSparseDiscreteBoundaryOperator']


# This imports dolfin at the same time as bempp if available to avoid delays
# at later imports of dolfin

try:
    import dolfin
except ImportError:
    have_dolfin = False  # pylint: disable=C0103
else:
    have_dolfin = True  # pylint: disable=C0103


# Check if config directory exists. If not create it.

def _check_create_init_dir():
    """Create the temporary dir if necessary."""
    from os.path import expanduser, join, isdir
    from os import mkdir

    home = expanduser("~")
    confpath = join(home, ".bempp")
    tmpath = join(confpath, "tmp")

    if not isdir(confpath):
        mkdir(confpath)
    if not isdir(tmpath):
        mkdir(tmpath)

    return confpath, tmpath


config_path, tmp_path = _check_create_init_dir()  # pylint: disable=C0103


# Get the path to Gmsh

def _gmsh_path():
    """Find Gmsh."""
    from .utils import which

    gmp = which("gmsh")
    if gmp is None:
        print("Could not find Gmsh. Interactive plotting and shapes module not available.")
    return gmp


gmsh_path = _gmsh_path()  # pylint: disable=C0103


# Define the global default options

from bempp.common import global_parameters as __global_parameters

global_parameters = __global_parameters()  # pylint: disable=C0103

# Now all the module imports

from bempp.grid import Grid, grid_from_sphere, grid_from_element_data, structured_grid
from bempp.grid import GridFactory
from bempp.space import function_space
from bempp.assembly import GridFunction
from bempp.assembly import InverseSparseDiscreteBoundaryOperator
from bempp import shapes
from bempp.file_interfaces import import_grid
from bempp import operators


def test():
    """ Runs BEM++ python unit tests """
    import unittest
    from os.path import dirname
    loader = unittest.TestLoader()
    suite = loader.discover(dirname(__file__))
    test_runner = unittest.TextTestRunner(verbosity=2)
    test_runner.run(suite)
