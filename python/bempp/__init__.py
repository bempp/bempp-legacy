""" Boundary Element Method package BEM++ """
from __future__ import print_function

__all__ = ['as_matrix',
           'grid_from_element_data',
           'structured_grid',
           'function_space',
           'GridFunction',
           'global_parameters',
           'GridFactory',
           'import_grid',
           'export'
           'operators',
           'shapes',
           'InverseSparseDiscreteBoundaryOperator',
           'ZeroBoundaryOperator']


# This imports dolfin at the same time as bempp if available to avoid delays
# at later imports of dolfin

try:
    import dolfin
except ImportError:
    HAVE_DOLFIN = False
else:
    HAVE_DOLFIN = True


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


CONFIG_PATH, TMP_PATH = _check_create_init_dir()


# Get the path to Gmsh

def _gmsh_path():
    """Find Gmsh."""
    from .utils import which

    gmp = which("gmsh")
    if gmp is None:
        print("Could not find Gmsh. Interactive plotting and shapes module not available.")
    return gmp


GMSH_PATH = _gmsh_path()


# Define the global default options

from bempp.common import global_parameters as __global_parameters


global_parameters = __global_parameters()

# Now all the module imports

from bempp.grid import grid_from_element_data, structured_grid
from bempp.grid import GridFactory
from bempp.space import function_space
from bempp.assembly import GridFunction
from bempp.assembly import InverseSparseDiscreteBoundaryOperator
from bempp.assembly import ZeroBoundaryOperator
from bempp.assembly import as_matrix
from bempp import shapes
from bempp.file_interfaces import import_grid
from bempp.file_interfaces import export
from bempp import operators


def test():
    """ Runs BEM++ python unit tests """
    import unittest
    from os.path import dirname
    loader = unittest.TestLoader()
    suite = loader.discover(dirname(__file__))
    test_runner = unittest.TextTestRunner(verbosity=2)
    test_runner.run(suite)
