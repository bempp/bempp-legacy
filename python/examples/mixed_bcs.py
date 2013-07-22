#!/usr/bin/env python

# This script solves a Laplace problem with mixed boundary conditions
# (part Dirichlet, part Neumann)

# Help Python find the bempp module
import sys
sys.path.append("..")

from bempp.lib import *
import numpy as np

# Boundary conditions (derived from an exact solution)

def evalDirichletData(point):
    x, y, z = point
    r = np.sqrt(x**2 + y**2 + z**2)
    return 2 * x * z / r**5 - y / r**3

def evalNeumannData(point):
    x, y, z = point
    r = np.sqrt(x**2 + y**2 + z**2)
    return -6 * x * z / r**6 + 2 * y / r**4

# Load mesh

grid = createGridFactory().importGmshGrid(
    "triangular", "../../examples/meshes/sphere-domains.msh")

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

# Initialize grid segments on which the Dirichlet and Neumann data are defined

segmentD = (GridSegment.closedDomain(grid, 1).
            union_(GridSegment.closedDomain(grid, 2)).
            union_(GridSegment.closedDomain(grid, 3)).
            union_(GridSegment.closedDomain(grid, 4)))
segmentN = segmentD.complement()

# Initialize spaces

pwiseConstants = createPiecewiseConstantScalarSpace(context, grid)
pwiseConstantsD = createPiecewiseConstantScalarSpace(context, grid, segmentD)
pwiseConstantsN = createPiecewiseConstantScalarSpace(context, grid, segmentN)
pwiseLinears = createPiecewiseLinearContinuousScalarSpace(context, grid)
pwiseLinearsD = createPiecewiseLinearContinuousScalarSpace(
    context, grid, segmentD)
pwiseLinearsN = createPiecewiseLinearContinuousScalarSpace(
    context, grid, segmentN)
# This space will be needed during the discretization of Dirichlet data
pwiseDLinearsD = createPiecewiseLinearDiscontinuousScalarSpace(
    context, grid, segmentD, strictlyOnSegment=True)

# Construct elementary operators

slpOpDD = createLaplace3dSingleLayerBoundaryOperator(
    context, pwiseConstantsD, pwiseLinearsD, pwiseConstantsD)
dlpOpDN = createLaplace3dDoubleLayerBoundaryOperator(
    context, pwiseLinearsN, pwiseLinearsD, pwiseConstantsD)
adlpOpND = createLaplace3dAdjointDoubleLayerBoundaryOperator(
    context, pwiseConstantsD, pwiseConstantsN, pwiseLinearsN)
hypOpNN = createLaplace3dHypersingularBoundaryOperator(
    context, pwiseLinearsN, pwiseConstantsN, pwiseLinearsN)

slpOpDN = createLaplace3dSingleLayerBoundaryOperator(
    context, pwiseConstantsN, pwiseLinearsD, pwiseConstantsD)
dlpOpDD = createLaplace3dDoubleLayerBoundaryOperator(
    context, pwiseLinearsD, pwiseLinearsD, pwiseConstantsD)
idOpDD = createIdentityOperator(
    context, pwiseLinearsD, pwiseLinearsD, pwiseConstantsD)
adlpOpNN = createLaplace3dAdjointDoubleLayerBoundaryOperator(
    context, pwiseConstantsN, pwiseConstantsN, pwiseLinearsN)
idOpNN = createIdentityOperator(
    context, pwiseConstantsN, pwiseConstantsN, pwiseLinearsN)
hypOpND = createLaplace3dHypersingularBoundaryOperator(
    context, pwiseLinearsD, pwiseConstantsN, pwiseLinearsN)

# Form the left- and right-hand-side operators

lhsOp = createBlockedBoundaryOperator(
    context,
    [[slpOpDD, -dlpOpDN],
     [adlpOpND, hypOpNN]])

# Construct the grid function representing the input data

dirichletData = createGridFunction(
    context, pwiseLinearsD, pwiseDLinearsD, evalDirichletData)
neumannData = createGridFunction(
    context, pwiseConstantsN, pwiseConstantsN, evalNeumannData)

# Construct the right-hand-side grid function

rhs = [(-0.5 * idOpDD + dlpOpDD) * dirichletData - slpOpDN * neumannData,
       -hypOpND * dirichletData + (-0.5 * idOpNN - adlpOpNN) * neumannData]

# Initialize the solver

solver = createDefaultIterativeSolver(lhsOp)
solver.initializeSolver(defaultGmresParameterList(1e-8))

# Solve the equation

solution = solver.solve(rhs)
print solution.solverMessage()

# Extract the solution in the form of grid functions

neumannSolution = solution.gridFunction(0)
dirichletSolution = solution.gridFunction(1)

# Combine known and unknown data

scatterPwiseConstantsD = createIdentityOperator(
    context, pwiseConstantsD, pwiseConstants, pwiseConstants)
scatterPwiseConstantsN = createIdentityOperator(
    context, pwiseConstantsN, pwiseConstants, pwiseConstants)
scatterPwiseLinearsD = createIdentityOperator(
    context, pwiseLinearsD, pwiseLinears, pwiseLinears)
scatterPwiseLinearsN = createIdentityOperator(
    context, pwiseLinearsN, pwiseLinears, pwiseLinears)

dirichlet = (scatterPwiseLinearsD * dirichletData +
             scatterPwiseLinearsN * dirichletSolution)
neumann = (scatterPwiseConstantsD * neumannSolution +
           scatterPwiseConstantsN * neumannData)

# Export data to VTK files

dirichlet.exportToVtk("vertex_data", "dirichlet_data", "dirichlet_data")
neumann.exportToVtk("cell_data", "neumann_data", "neumann_data")

# Compare the numerical and analytical solution

evalOptions = createEvaluationOptions()
absError, relError = estimateL2Error(dirichlet, evalDirichletData,
                                     quadStrategy, evalOptions)
print "Relative L^2 error (Dirichlet):", relError

evalOptions = createEvaluationOptions()
absError, relError = estimateL2Error(neumann, evalNeumannData,
                                     quadStrategy, evalOptions)
print "Relative L^2 error (Neumann):", relError
