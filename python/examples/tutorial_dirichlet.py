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

# Create quadrature strategy

accuracyOptions = createAccuracyOptions()
# Increase by 2 the order of quadrature rule used to approximate
# integrals of regular functions on pairs on elements
accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2)
# Increase by 2 the order of quadrature rule used to approximate
# integrals of regular functions on single elements
accuracyOptions.singleRegular.setRelativeQuadratureOrder(2)
quadStrategy = createNumericalQuadratureStrategy(
    "float64", "float64", accuracyOptions)

# Create assembly context

assemblyOptions = createAssemblyOptions()
assemblyOptions.switchToAcaMode(createAcaOptions())
context = createContext(quadStrategy, assemblyOptions)

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

# Extract the solution in the form of a grid function and export it
# in VTK format

solFun = solution.gridFunction()
solFun.exportToVtk("cell_data", "neumann_data", "solution")

# Compare the numerical and analytical solution

evalOptions = createEvaluationOptions()
absError, relError = estimateL2Error(solFun, evalExactNeumannData,
                                     quadStrategy, evalOptions)
print "Relative L^2 error:", relError

# Prepare to evaluate the solution on an annulus outside the sphere

# Create potential operators

slPotOp = createLaplace3dSingleLayerPotentialOperator(context)
dlPotOp = createLaplace3dDoubleLayerPotentialOperator(context)

# Define points at which the solution should be evaluated

rCount = 51;
thetaCount = 361;
r, theta, z = np.mgrid[1:2:rCount*1j, 0:2*np.pi:thetaCount*1j, 0:0:1j]
x = r * np.cos(theta)
y = r * np.sin(theta)
# put the x, y and z coordinates in successive rows of a matrix
evaluationPoints = np.vstack((x.ravel(), y.ravel(), z.ravel()))

# Use the Green's representation formula to evaluate the solution

evaluationOptions = createEvaluationOptions()
field = (-slPotOp.evaluateAtPoints(solFun, evaluationPoints,
                                   evaluationOptions) +
          dlPotOp.evaluateAtPoints(dirichletData, evaluationPoints,
                                   evaluationOptions))

# Plot data

from bempp import visualization2 as vis
annulus = vis.tvtkStructuredGridData(
    evaluationPoints, field, (rCount, thetaCount))
sphere = vis.tvtkGridFunction(dirichletData)
vis.plotScalarData(tvtkGridFunctions=sphere,
                   tvtkStructuredGridData=annulus)
