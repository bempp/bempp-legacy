#!/usr/bin/env python

# This script solves the Maxwell equations in the region exterior to a bounded
# object, with Dirichlet boundary conditions given by the exact solution
# (satisfying the Silver-Mueller radiation conditions)
#
#     \vec u(\vec x) = h_1^{(1)}(k r) \hat phi,
#
# where (r, theta, phi) are the radial, zenith angle and azimuthal spherical
# coordinates in the system anchored at the point (0.1, 0.1, 0.1), h_1^{(1)}(r)
# is the spherical Hankel function of the first kind and order 1 and \hat phi is
# the unit vector oriented along d(\vec x)/d\phi.

# Help Python find the bempp module
import sys
sys.path.append("..")

from bempp.lib import *
import numpy as np
from cmath import exp
from math import sqrt

# Parameters

k = 2
source = 0.1

# Boundary conditions

def evalDirichletData(point, normal):
    x, y, z = point - source
    r = sqrt(x**2 + y**2 + z**2)
    kr = k * r
    h1kr = (-1j - kr) * exp(1j * kr) / (kr * kr)
    scale = h1kr / r
    field = [-y * scale, x * scale, 0.]
    return np.cross(field, normal)

# Exact solution

def evalExactNeumannData(point, normal):
    x, y, z = point - source
    r = sqrt(x**2 + y**2 + z**2)
    kr = k * r
    h1kr = (-1j - kr) * exp(1j * kr) / (kr * kr)
    h1kr_deriv = ((1. + 1j - 1j * kr) * (1. + 1j + kr) *
                  exp(1j * kr) / (kr * kr * r))
    xy_factor = (h1kr - r * h1kr_deriv) / (r * r * r)
    curl = [x * z * xy_factor,
            y * z * xy_factor,
            ((x*x + y*y + 2*z*z) * h1kr + r * (x*x + y*y) * h1kr_deriv) /
            (r * r * r)]
    return np.cross(curl, normal) / (1j * k)

def evalExactSolution(point):
    x, y, z = point - source
    r = sqrt(x**2 + y**2 + z**2)
    kr = k * r
    h1kr = (-1j - kr) * exp(1j * kr) / (kr * kr)
    scale = h1kr / r
    return np.array([-y * scale, x * scale, 0.])

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
    "float64", "complex128", accuracyOptions)

# Create assembly context

assemblyOptions = createAssemblyOptions()
assemblyOptions.switchToAcaMode(createAcaOptions())
context = createContext(quadStrategy, assemblyOptions)

# Initialize spaces

space = createRaviartThomas0VectorSpace(context, grid)

# Construct elementary operators

slpOp = createMaxwell3dSingleLayerBoundaryOperator(
    context, space, space, space, k, "SLP")
dlpOp = createMaxwell3dDoubleLayerBoundaryOperator(
    context, space, space, space, k, "DLP")
idOp = createMaxwell3dIdentityOperator(
    context, space, space, space, "Id")

# Form the left- and right-hand-side operators

lhsOp = slpOp
rhsOp = -(0.5 * idOp + dlpOp)

# Construct the grid function representing the (input) Dirichlet data

dirichletData = createGridFunction(
    context, space, space, evalDirichletData, True)

# Construct the right-hand-side grid function

rhs = rhsOp * dirichletData

# Initialize the solver

solver = createDefaultIterativeSolver(lhsOp)
precTol = 1e-2
invLhsOp = acaOperatorApproximateLuInverse(
    lhsOp.weakForm().asDiscreteAcaBoundaryOperator(), precTol)
prec = discreteOperatorToPreconditioner(invLhsOp)
solver.initializeSolver(defaultGmresParameterList(1e-8), prec)

# Solve the equation

solution = solver.solve(rhs)
print solution.solverMessage()

# Extract the solution in the form of a grid function and export it
# in VTK format

neumannData = solution.gridFunction()
neumannData.exportToVtk("vertex_data", "neumann_data", "calculated_neumann_data_vertex")

# Export the exact solution (for comparison)
exactNeumannData = createGridFunction(
    context, space, space, evalExactNeumannData, True)
exactNeumannData.exportToVtk("vertex_data", "neumann_data", "exact_neumann_data_vertex")

# Compare the numerical and analytical solution

evalOptions = createEvaluationOptions()
absError, relError = estimateL2Error(neumannData, evalExactNeumannData,
                                     quadStrategy, evalOptions, True)
print "Relative L^2 error:", relError

# Prepare to evaluate the solution on an annulus outside the sphere

# Create potential operators

slPotOp = createMaxwell3dSingleLayerPotentialOperator(context, k)
dlPotOp = createMaxwell3dDoubleLayerPotentialOperator(context, k)

# TODO: make the code below more interesting and add some plotting; the current
# version just demonstrates that the potential operators are defined correctly

# Define points at which the solution should be evaluated

x = np.array([3])
y = np.array([2])
z = np.array([1])
# put the x, y and z coordinates in successive rows of a matrix
evaluationPoints = np.vstack((x.ravel(), y.ravel(), z.ravel()))

# Use the Green's representation formula to evaluate the solution

evaluationOptions = createEvaluationOptions()
field = (-slPotOp.evaluateAtPoints(neumannData, evaluationPoints,
                                   evaluationOptions)
         -dlPotOp.evaluateAtPoints(dirichletData, evaluationPoints,
                                   evaluationOptions))
exactField = evalExactSolution(np.array([3, 2, 1]))
print "field:", field.ravel()
print "exact field:", exactField
