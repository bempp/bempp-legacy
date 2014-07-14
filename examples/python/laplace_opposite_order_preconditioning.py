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
    "triangular", "../../meshes/sphere-h-0.2.msh")

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

# The space of constant functions is defined on the dual mesh. This gives
# a stable pairing between the piecewise linear and the constant functions.
pwiseConstants = createPiecewiseConstantDualGridScalarSpace(context, grid)
pwiseLinears = createPiecewiseLinearContinuousScalarSpace(context, grid)

# Construct elementary operators

slpOp = createLaplace3dSingleLayerBoundaryOperator(
    context, pwiseConstants, pwiseLinears, pwiseConstants)
dlpOp = createLaplace3dDoubleLayerBoundaryOperator(
    context, pwiseLinears, pwiseLinears, pwiseConstants)
idOp = createIdentityOperator(
    context, pwiseLinears, pwiseLinears, pwiseConstants)
hypOp = createLaplace3dHypersingularBoundaryOperator(
    context, pwiseLinears, pwiseConstants, pwiseLinears)

# Form the left- and right-hand-side operators
# Note the multiplication with the hypersingular operator to precondition
# the problem.

# No opposite order preconditioning
lhsOp1 = slpOp
rhsOp1 = -0.5 * idOp + dlpOp


# Opposite order preconditioning
lhsOp2 = hypOp*(slpOp)
rhsOp2 = hypOp*(-0.5 * idOp + dlpOp)

# Construct the grid function representing the (input) Dirichlet data

dirichletData = createGridFunction(
    context, pwiseLinears, pwiseLinears, evalDirichletData)

# Construct the right-hand-side grid function

rhs1 = rhsOp1 * dirichletData
rhs2 = rhsOp2 * dirichletData

# Initialize the solvers

# Standard solver
solver1 = createDefaultIterativeSolver(lhsOp1)
solver1.initializeSolver(defaultGmresParameterList(1e-5))

# Solve with mass matrix preconditioning
solver2= createDefaultIterativeSolver(lhsOp1,"test_convergence_in_range")
solver2.initializeSolver(defaultGmresParameterList(1e-5))

# Solve the full opposite order formulation
solver3 = createDefaultIterativeSolver(lhsOp2,"test_convergence_in_range")
solver3.initializeSolver(defaultGmresParameterList(1e-5))



# Solve the equation

solution1 = solver1.solve(rhs1)
solution2 = solver2.solve(rhs1)
solution3 = solver3.solve(rhs2)


print "Results using GMRES"
print 

print "Number of iterations for unpreconditioned solve: "+str(solution1.iterationCount())
print "Number of iterations for mass matrix preconditioned solve: "+str(solution2.iterationCount())
print "Number of iterations for opposite order formulation: "+str(solution3.iterationCount())


