""" Boundary Element Method package BEM++ """
from __future__ import print_function

# This imports dolfin at the same time as bempp if available to avoid delays
# at later imports of dolfin

try:
    import dolfin as _
except ImportError:
    HAVE_DOLFIN = False
else:
    HAVE_DOLFIN = True

# Initialize logger

from bempp.utils.logging import _init_logger

LOGGER = _init_logger()

# Check if config directory exists. If not create it.

def _check_create_init_dir():
    """Create the temporary dir if necessary."""
    from os.path import expanduser, join, isdir
    from os import mkdir

    home = expanduser("~")
    config_path = join(home,".bempp")

    try:
        if not isdir(config_path):
            mkdir(config_path)
    except OSError: # Read only file system try a tmp dir
        import tempfile
        import warnings
        warnings.warn("Could not create BEM++ config dir."
            "Falling back to a temorary dir."
            "Your config will not be stored")
        config_path = tempfile.mkdtemp()

    tmp_path = join(config_path,"tmp")
    if not isdir(tmp_path):
        mkdir(tmp_path)

    return config_path, tmp_path


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
from bempp.assembly import assemble_dense_block
from bempp import shapes
from bempp.file_interfaces import import_grid
from bempp.file_interfaces import export
from bempp import operators

from bempp.utils.logging import DEBUG, INFO, WARNING, ERROR, CRITICAL
from bempp.utils.logging import enable_console_logging
from bempp.utils.logging import enable_file_logging
from bempp.utils.logging import set_logging_level

ALL = -1 # Useful global identifier

def test():
    """ Runs BEM++ python unit tests """
    import unittest
    from os.path import dirname
    loader = unittest.TestLoader()
    suite = loader.discover(dirname(__file__))
    test_runner = unittest.TextTestRunner(verbosity=2)
    test_runner.run(suite)
