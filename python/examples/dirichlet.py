#!/usr/bin/env python

# This script imports a grid in the Gmsh format and exports the z coordinate 
# of each vertex to the file "output.vtu" in the VTK format.

# Help Python find the bempp module
import sys
sys.path.append("..")

from bempp.lib import *
import numpy as np


def evalDirichletData(point):
    return -1

print "Importing grid..."
grid = createGridFactory().importGmshGrid("triangular",
                                        "../../examples/meshes/sphere-152.msh")

factory = createNumericalQuadratureStrategy("float64", "float64")
options = createAssemblyOptions()
options.switchToAca(createAcaOptions())

context = createContext(factory, options)

pwiseConstants = createPiecewiseConstantScalarSpace(context, grid)
pwiseLinears = createPiecewiseLinearContinuousScalarSpace(context, grid)
pwiseConstants.assignDofs()
pwiseLinears.assignDofs()

slpOp = createLaplace3dSingleLayerBoundaryOperator(
    context, pwiseConstants, pwiseLinears, pwiseConstants)
dlpOp = createLaplace3dDoubleLayerBoundaryOperator(
    context, pwiseLinears, pwiseLinears, pwiseConstants)
idOp = createIdentityOperator(
    context, pwiseLinears, pwiseLinears, pwiseConstants)

lhsOp = slpOp
rhsOp = -0.5 * idOp + dlpOp

print "Evaluating Dirichlet data..."
dirichletData = gridFunctionFromSurfaceNormalIndependentFunction(
    context, pwiseLinears, pwiseLinears, evalDirichletData)

rhs = rhsOp * dirichletData

solver = createDefaultIterativeSolver(lhsOp)
params = defaultGmresParameterList(1e-5)
solver.initializeSolver(params)
solution = solver.solve(rhs)
neumannData = solution.gridFunction()

neumannData.exportToVtk("vertex_data","neumann_data", "neumann_data_vertex")
neumannData.exportToVtk("cell_data","neumann_data", "neumann_data_cell")
