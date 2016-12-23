"""This module contains basic classes for the assembly of integral operators."""

from bempp.api.assembly.discrete_boundary_operator import \
    GeneralNonlocalDiscreteBoundaryOperator
from bempp.api.assembly.discrete_boundary_operator import \
    DenseDiscreteBoundaryOperator
from bempp.api.assembly.discrete_boundary_operator import \
    SparseDiscreteBoundaryOperator
from bempp.api.assembly.discrete_boundary_operator import \
    InverseSparseDiscreteBoundaryOperator
from bempp.api.assembly.discrete_boundary_operator import \
    ZeroDiscreteBoundaryOperator
from bempp.api.assembly.discrete_boundary_operator import \
    DiscreteRankOneOperator
from bempp.api.assembly.discrete_boundary_operator import as_matrix
from bempp.api.assembly.boundary_operator import BoundaryOperator
from bempp.api.assembly.boundary_operator import LocalBoundaryOperator
from bempp.api.assembly.boundary_operator import ElementaryBoundaryOperator
from bempp.api.assembly.boundary_operator import ZeroBoundaryOperator
from bempp.api.assembly.boundary_operator import RankOneBoundaryOperator
from bempp.api.assembly.blocked_operator import BlockedOperator
from bempp.api.assembly.blocked_operator import BlockedDiscreteOperator
from bempp.api.assembly.grid_function import GridFunction
from bempp.api.assembly.assembler import assemble_dense_block
from bempp.api.assembly.potential_operator import PotentialOperator

from bempp.api.assembly import functors
