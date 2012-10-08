#!/usr/bin/env python

# Help Python find the bempp module
import sys
sys.path.append("..")

from bempp.lib import *
import numpy as np

# Boundary conditions

def evalDirichletData(point):
    x, y, z = point
    r = np.sqrt(x**2 + y**2 + z**2)
    return 2 * x * z / r**5 - y / r**3

# Exact solution

def evalExactNeumannData(point):
    x, y, z = point
    r = np.sqrt(x**2 + y**2 + z**2)
    return -6 * x * z / r**6 + 2 * y / r**4

# Load mesh

grid = createGridFactory().importGmshGrid(
    "triangular", "../../examples/meshes/sphere-h-0.2.msh")

# Create assembly context

quadStrategy = createNumericalQuadratureStrategy("float64", "float64")
options = createAssemblyOptions()
options.switchToAcaMode(createAcaOptions())
context = createContext(quadStrategy, options)

# Initialize spaces

pwiseConstants = createPiecewiseConstantScalarSpace(context, grid)
pwiseLinears = createPiecewiseLinearContinuousScalarSpace(context, grid)

# Construct elementary operators

slpOp = createLaplace3dSingleLayerBoundaryOperator(
    context, pwiseConstants, pwiseLinears, pwiseConstants)
dlpOp = createLaplace3dDoubleLayerBoundaryOperator(
    context, pwiseLinears, pwiseLinears, pwiseConstants)
idOp = createIdentityOperator(
    context, pwiseLinears, pwiseLinears, pwiseConstants)

# Form the left- and right-hand-side operators

lhsOp = slpOp
rhsOp = -0.5 * idOp + dlpOp

# Construct the grid function representing the (input) Dirichlet data

dirichletData = createGridFunction(
    context, pwiseLinears, pwiseLinears, evalDirichletData)

# Construct the right-hand-side grid function

rhs = rhsOp * dirichletData

# Initialize the solver

solver = createDefaultIterativeSolver(lhsOp)
solver.initializeSolver(defaultGmresParameterList(1e-5))

# Solve the equation

solution = solver.solve(rhs)
print solution.solverMessage()

# Extract the solution in the form of a grid function and export it in VTK format

solFun = solution.gridFunction()
solFun.exportToVtk("cell_data", "neumann_data", "solution")

# Compare the numerical and analytical solution

exactSolFun = createGridFunction(
    context, pwiseConstants, pwiseConstants, evalExactNeumannData)
diff = solFun - exactSolFun

relError = diff.L2Norm() / exactSolFun.L2Norm()
print "Relative L^2 error:", relError

