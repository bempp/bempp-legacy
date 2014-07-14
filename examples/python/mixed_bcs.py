#!/usr/bin/env python

# This script solves a Laplace problem with mixed boundary conditions
# (part Dirichlet, part Neumann)

# Help Python find the bempp module
import sys
sys.path.append("..")

from bempp.lib import *
import numpy as np

# Load mesh. The chosen mesh is a sphere divided into 8 domains (with indices 1
# to 8).

grid = createGridFactory().importGmshGrid(
    "triangular", "../../meshes/sphere-domains.msh")

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

# Dirichlet boundary conditions will be imposed on domains 1 to 4, including
# their boundary
segmentD = (GridSegment.closedDomain(grid, 1)
            .union_(GridSegment.closedDomain(grid, 2))
            .union_(GridSegment.closedDomain(grid, 3))
            .union_(GridSegment.closedDomain(grid, 4)))
# Neumann boundary conditions will be imposed on the rest of the grid
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
# This space will be needed during the discretization of Dirichlet data.
# strictlyOnSegment=True means that the basis functions are truncated to the
# elements that belong to the segment. So, for example, a function associated
# with a vertex lying at the boundary of the segment will be set to zero on all
# elements not belonging to the segment, hence in fact discontinuous on the grid
# as a whole (but not on the segment itself).
pwiseDLinearsD = createPiecewiseLinearContinuousScalarSpace(
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

# Note that we pass pwiseDLinearsD as the dual space when discretizing the
# Dirichlet data. You can see for yourself that if pwiseLinearsD is passed
# there, the approximation quality of the Dirichlet data near the boundary of
# the segment is degraded and as a consequence the solution has a much larger
# L^2 error.
dirichletData = createGridFunction(
    context, pwiseLinearsD, pwiseDLinearsD, evalDirichletData)
neumannData = createGridFunction(
    context, pwiseConstantsN, pwiseConstantsN, evalNeumannData)

# Construct the right-hand-side grid functions

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

# Combine imposed and calculated parts of traces into single grid functions

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

