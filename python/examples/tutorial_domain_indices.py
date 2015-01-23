#!/usr/bin/env python

# This script solves the Laplace equation in the space outside
# two unit spheres, one put at potential +1, the other -1.
# It demonstrates the construction of grid functions dependent
# on domain indices from the Gmsh file.

# Help Python find the bempp module
import sys
sys.path.append("..")

from bempp.lib import *
import numpy as np

# Boundary conditions

def evalDirichletData(point, domain_index):
    if domain_index == 1:
        return 1.
    elif domain_index == 2:
        return -1.
    else:
        raise ValueError("Unexpected domain index")

# Load mesh

grid = createGridFactory().importGmshGrid(
    "triangular", "../../examples/meshes/two-spheres.msh")

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
    context, pwiseLinears, pwiseLinears, evalDirichletData,
    domainIndexDependent=True)

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

# Prepare to evaluate the solution on a fragment of the xz plane

# Create potential operators

slPotOp = createLaplace3dSingleLayerPotentialOperator(context)
dlPotOp = createLaplace3dDoubleLayerPotentialOperator(context)

# Define points at which the solution should be evaluated

xCount = 101;
yCount = 61;
x, y, z = np.mgrid[-5:5:xCount*1j, -3:3:yCount*1j, 0:0:1j]
# put the x, y and z coordinates in successive rows of a matrix
evaluationPoints = np.vstack((x.ravel(), y.ravel(), z.ravel()))

# Select the points located outside the two spheres

outside = np.logical_not(areInside(grid, evaluationPoints))

# Use the Green's representation formula to evaluate the solution

evaluationOptions = createEvaluationOptions()
fieldOutside = (-slPotOp.evaluateAtPoints(solFun,
                                          evaluationPoints[:,outside],
                                          evaluationOptions) +
                 dlPotOp.evaluateAtPoints(dirichletData,
                                          evaluationPoints[:,outside],
                                          evaluationOptions))

field = np.zeros(xCount * yCount)
np.place(field, outside, fieldOutside.ravel())

# Plot data

from bempp import visualization2 as vis
visField = vis.tvtkStructuredGridData(
    evaluationPoints, field, (xCount, yCount))
visDirichlet = vis.tvtkGridFunction(dirichletData)
vis.plotScalarData(tvtkGridFunctions=visDirichlet,
                   tvtkStructuredGridData=visField)
