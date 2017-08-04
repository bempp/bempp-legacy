""" Boundary Element Method package BEM++ """
from __future__ import print_function
import os

from mpi4py import MPI


#pylint: disable=bare-except
#pylint: disable=wrong-import-position
#pylint: disable=invalid-name
#pylint: disable=ungrouped-imports

mpi_comm = MPI.COMM_WORLD #pylint: disable=no-member
mpi_rank = mpi_comm.Get_rank() #pylint: disable=no-member
mpi_size = mpi_comm.Get_size() #pylint: disable=no-member
mpi_name = MPI.Get_processor_name() #pylint: disable=no-member

PLOT_BACKEND = 'gmsh'

# import the version string
from bempp import config as _config
__version__ = _config.version


# Initialize logger

from bempp.api.utils.logging import _init_logger
from bempp.api.utils.logging import log

_console_logging_handler = None
_LOGGER = _init_logger()


try:
    if os.environ['BEMPP_CONSOLE_LOGGING'] == '1':
        from bempp.api.utils.logging import enable_console_logging
        enable_console_logging()
except:
    pass

# Check for FEniCS

try:
    #pylint: disable=import-error
    import dolfin as _
except:
    HAVE_DOLFIN = False
    log(
            "Dolfin could not be imported." +
            "FEM/BEM coupling with FEniCS not available.")
else:
    HAVE_DOLFIN = True
    log("Found Dolfin. FEM/BEM coupling with FEniCS enabled.")


# Check if config directory exists. If not create it.

def _check_create_init_dir():
    """Create the temporary dir if necessary."""
    from os.path import expanduser, join, isdir
    from os import mkdir
    import tempfile

    home = expanduser("~")
    config_path = join(home, ".bempp")

    try:
        if not isdir(config_path):
            mkdir(config_path)
    except OSError:  # Read only file system try a tmp dir
        import warnings
        warnings.warn("Could not create BEM++ config dir."
                "Falling back to a temorary dir."
                "Your config will not be stored")
        config_path = tempfile.mkdtemp()

    tmp_path = tempfile.mkdtemp()

    return config_path, tmp_path


CONFIG_PATH, TMP_PATH = _check_create_init_dir()

# Get the path to Gmsh


def _gmsh_path():
    """Find Gmsh."""
    from bempp.api.utils import which

    gmp = which("gmsh")
    if gmp is None:
        print(
                "Could not find Gmsh." +
                "Interactive plotting and shapes module not available.")
    return gmp


GMSH_PATH = _gmsh_path()


# Define the global default options

from bempp.api.common import global_parameters as __global_parameters
global_parameters = __global_parameters()

# Now all the module imports

from bempp.api.grid import grid_from_element_data, structured_grid
from bempp.api.grid import GridFactory
from bempp.api.space import function_space
from bempp.api.space import project_operator
from bempp.api.assembly import GridFunction
from bempp.api.assembly import InverseSparseDiscreteBoundaryOperator
from bempp.api.assembly import ZeroBoundaryOperator
from bempp.api.assembly import RankOneBoundaryOperator
from bempp.api.assembly import as_matrix
from bempp.api.assembly import assemble_dense_block
from bempp.api.assembly import BlockedOperator
from bempp.api.assembly import BlockedDiscreteOperator
from bempp.api.assembly import functors
from bempp.api import shapes
from bempp.api.file_interfaces import import_grid
from bempp.api.file_interfaces import export
from bempp.api.file_interfaces import three_planes_view
from bempp.api.external.viewers import set_gmsh_viewer
from bempp.api.external.viewers import set_ipython_notebook_viewer
from bempp.api import operators
from bempp.api import linalg
from bempp.api import hmat
from bempp.api import _fmm
from bempp.api.hmat import hmatrix_interface

from bempp.api.utils.logging import enable_console_logging
from bempp.api.utils.logging import flush_log
from bempp.api.utils.logging import enable_file_logging
from bempp.api.utils.logging import set_logging_level
from bempp.api.utils.logging import timeit
from bempp.api.utils.logging import Timer

ALL = -1  # Useful global identifier


def test():
    """ Runs BEM++ python unit tests """
    import unittest
    from os.path import dirname
    loader = unittest.TestLoader()
    suite = loader.discover(dirname(__file__))
    test_runner = unittest.TextTestRunner(verbosity=2)
    result = test_runner.run(suite)
    return result.wasSuccessful()
