"""This module contains basic classes for the assembly of integral operators."""

from .discrete_boundary_operator import GeneralNonlocalDiscreteBoundaryOperator
from .discrete_boundary_operator import DenseDiscreteBoundaryOperator
from .discrete_boundary_operator import SparseDiscreteBoundaryOperator
from .discrete_boundary_operator import InverseSparseDiscreteBoundaryOperator
from .boundary_operator import LocalBoundaryOperator
from .boundary_operator import ElementaryBoundaryOperator
from .boundary_operator import ZeroBoundaryOperator
from .grid_function import GridFunction



