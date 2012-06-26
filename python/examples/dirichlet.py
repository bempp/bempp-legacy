#!/usr/bin/env python

# This script imports a grid in the Gmsh format and exports the z coordinate 
# of each vertex to the file "output.vtu" in the VTK format.

# Help Python find the bempp module
import sys
sys.path.append("..")

import bempp
import numpy as np


def eval_dirichlet_data(point):
    return -1

print "Importing grid..."
grid = bempp.GridFactory.importGmshGrid("triangular",
                                        "../../examples/meshes/double-sphere-5162.msh")

pwiseConstants = bempp.piecewiseConstantScalarSpace(grid)
pwiseLinears = bempp.piecewiseLinearContinuousScalarSpace(grid)
pwiseConstants.assignDofs()
pwiseLinears.assignDofs()

slp = bempp.laplace3dSingleLayerPotential(pwiseConstants, pwiseConstants)
dlp = bempp.laplace3dDoubleLayerPotential(pwiseConstants, pwiseLinears)
id = bempp.identityOperator(pwiseConstants, pwiseLinears)

lhs_op = slp
rhs_op = -0.5 * id + dlp

factory = bempp.standardLocalAssemblerFactoryForOperatorsOnSurfaces()
options = bempp.AssemblyOptions()
options.switchToAca(bempp.AcaOptions())
print "Assembling LHS operator..."
lhs_op.assembleWeakForm(factory, options)
print "Assembling RHS operator..."
rhs_op.assembleWeakForm(factory, options)

print "Evaluating Dirichlet data..."
dirichlet_data_functor = bempp.surfaceNormalIndependentFunctor(
    eval_dirichlet_data)
dirichlet_data = bempp.gridFunctionFromSurfaceNormalIndependentFunctor(
    pwiseLinears, dirichlet_data_functor, factory, options)

rhs = rhs_op * dirichlet_data

solver = bempp.defaultIterativeSolver(lhs_op, rhs)
params = bempp.defaultGmresParameterList(1e-5)
solver.initializeSolver(params)
solver.solve()

neumann_data = solver.getResult()

neumann_data.exportToVtk(bempp.VtkWriter.VERTEX_DATA,
    "neumann_data", "neumann_data_vertex")
neumann_data.exportToVtk(bempp.VtkWriter.CELL_DATA,
    "neumann_data", "neumann_data_cell")

neumann_data.plot(mode=bempp.VtkWriter.CELL_DATA)
