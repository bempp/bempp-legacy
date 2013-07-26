#!/usr/bin/env python

# This script solves a Laplace problem with mixed boundary conditions
# (part Dirichlet, part Neumann) using higher-order basis functions

# Help Python find the bempp module
import sys
sys.path.append("..")

from bempp.lib import *
import numpy as np

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

orderN = 1
orderD = 2
neumannSpace = createPiecewisePolynomialDiscontinuousScalarSpace(
    context, grid, orderN,
    requireElementOnSegment=True, requireReferencePointOnSegment=False)
neumannSpaceD = createPiecewisePolynomialDiscontinuousScalarSpace(
    context, grid, orderN, segmentD,
    requireElementOnSegment=True, requireReferencePointOnSegment=False)
neumannSpaceN = createPiecewisePolynomialDiscontinuousScalarSpace(
    context, grid, orderN, segmentN,
    requireElementOnSegment=True, requireReferencePointOnSegment=False)
dirichletSpace = createPiecewisePolynomialContinuousScalarSpace(
    context, grid, orderD)
dirichletSpaceD = createPiecewisePolynomialContinuousScalarSpace(
    context, grid, orderD, segmentD)
dirichletSpaceN = createPiecewisePolynomialContinuousScalarSpace(
    context, grid, orderD, segmentN)
dualDirichletSpaceD = createPiecewisePolynomialContinuousScalarSpace(
    context, grid, orderD, segmentD, strictlyOnSegment=True)

# Construct elementary operators

slpOpDD = createLaplace3dSingleLayerBoundaryOperator(
    context, neumannSpaceD, dirichletSpaceD, neumannSpaceD)
dlpOpDN = createLaplace3dDoubleLayerBoundaryOperator(
    context, dirichletSpaceN, dirichletSpaceD, neumannSpaceD)
adlpOpND = createLaplace3dAdjointDoubleLayerBoundaryOperator(
    context, neumannSpaceD, neumannSpaceN, dirichletSpaceN)
hypOpNN = createLaplace3dHypersingularBoundaryOperator(
    context, dirichletSpaceN, neumannSpaceN, dirichletSpaceN)

slpOpDN = createLaplace3dSingleLayerBoundaryOperator(
    context, neumannSpaceN, dirichletSpaceD, neumannSpaceD)
dlpOpDD = createLaplace3dDoubleLayerBoundaryOperator(
    context, dirichletSpaceD, dirichletSpaceD, neumannSpaceD)
idOpDD = createIdentityOperator(
    context, dirichletSpaceD, dirichletSpaceD, neumannSpaceD)
adlpOpNN = createLaplace3dAdjointDoubleLayerBoundaryOperator(
    context, neumannSpaceN, neumannSpaceN, dirichletSpaceN)
idOpNN = createIdentityOperator(
    context, neumannSpaceN, neumannSpaceN, dirichletSpaceN)
hypOpND = createLaplace3dHypersingularBoundaryOperator(
    context, dirichletSpaceD, neumannSpaceN, dirichletSpaceN)

# Form the left-hand-side operator

lhsOp = createBlockedBoundaryOperator(
    context,
    [[slpOpDD, -dlpOpDN],
     [adlpOpND, hypOpNN]])

# Construct the grid functions representing the known parts of the Dirichlet and
# Neumann traces. They are derived from an exact solution of the Laplace
# equation

def evalDirichletData(point):
    x, y, z = point
    r = np.sqrt(x**2 + y**2 + z**2)
    return 2 * x * z / r**5 - y / r**3

def evalNeumannData(point):
    x, y, z = point
    r = np.sqrt(x**2 + y**2 + z**2)
    return -6 * x * z / r**6 + 2 * y / r**4

dirichletData = createGridFunction(
    context, dirichletSpaceD, dualDirichletSpaceD, evalDirichletData)
neumannData = createGridFunction(
    context, neumannSpaceN, neumannSpaceN, evalNeumannData)

# Construct the right-hand-side grid functions

rhs = [(-0.5 * idOpDD + dlpOpDD) * dirichletData - slpOpDN * neumannData,
       -hypOpND * dirichletData + (-0.5 * idOpNN - adlpOpNN) * neumannData]

# Initialize the solver

solver = createDefaultIterativeSolver(lhsOp)
solver.initializeSolver(defaultGmresParameterList(1e-8, 10000))

# Solve the equation

solution = solver.solve(rhs)
print solution.solverMessage()

# Extract the solution in the form of grid functions

neumannSolution = solution.gridFunction(0)
dirichletSolution = solution.gridFunction(1)

# Combine imposed and calculated parts of traces into single grid functions

scatterPwiseConstantsD = createIdentityOperator(
    context, neumannSpaceD, neumannSpace, neumannSpace)
scatterPwiseConstantsN = createIdentityOperator(
    context, neumannSpaceN, neumannSpace, neumannSpace)
scatterPwiseLinearsD = createIdentityOperator(
    context, dirichletSpaceD, dirichletSpace, dirichletSpace)
scatterPwiseLinearsN = createIdentityOperator(
    context, dirichletSpaceN, dirichletSpace, dirichletSpace)

dirichlet = (scatterPwiseLinearsD * dirichletData +
             scatterPwiseLinearsN * dirichletSolution)
neumann = (scatterPwiseConstantsD * neumannSolution +
           scatterPwiseConstantsN * neumannData)

# Export data to Gmsh files

dirichlet.exportToGmsh("dirichlet_data", "dirichlet_data.msh")
neumann.exportToGmsh("neumann_data", "neumann_data.msh")

# Compare the numerical and analytical solution

evalOptions = createEvaluationOptions()
absError, relError = estimateL2Error(dirichlet, evalDirichletData,
                                     quadStrategy, evalOptions)
print "Relative L^2 error (Dirichlet):", relError

evalOptions = createEvaluationOptions()
absError, relError = estimateL2Error(neumann, evalNeumannData,
                                     quadStrategy, evalOptions)
print "Relative L^2 error (Neumann):", relError

