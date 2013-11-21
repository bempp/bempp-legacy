# Scattering from a sphere using a combined direct formulation

# Help Python find the bempp module
import sys
sys.path.append("..")

# Import numpy and library modules

import numpy as np
from bempp import lib
# from bempp import visualization as vis

# Define the wavenumber

k = 5

def evalDirichletData(point):
    x, y, z = point
    r = np.sqrt(x**2 + y**2 + z**2)
    return x# 2 * x * z / r**5 - y / r**3

# Set accuracy options. For regular integrals on pairs on elements,
# use quadrature of 2 orders higher than default, and for singular
# integrals, 1 order higher than default.

accuracyOptions = lib.createAccuracyOptions()
accuracyOptions.doubleRegular.setRelativeQuadratureOrder(2)
accuracyOptions.doubleSingular.setRelativeQuadratureOrder(1)

# The default strategy for numerical quadrature.
# "float64" is the basis function type and "complex128" is the result type".

quadStrategy = lib.createNumericalQuadratureStrategy(
    "float64", "complex128", accuracyOptions)

# Use ACA to accelerate the assembly

options = lib.createAssemblyOptions()
#print dir(options)

if 0: # Use ACA to accelerate the assembly
	#options.switchToAca(lib.createAcaOptions())
	options.switchToDense()
else: # Use FMM to accelerate the assembly
	fmmOptions = lib.createFmmOptions()
	fmmOptions.levels = 4
	fmmOptions.expansionOrder = 3 # number of terms in addition theorem in leaves
	#fmmOptions.expansionOrderMax = 6 # max number of terms in addition theorem
	options.switchToFmmMode(fmmOptions)

# The context object combines these settings

context = lib.createContext(quadStrategy, options)

# Load a grid approximating a unit sphere

grid_factory = lib.createGridFactory()
grid = grid_factory.importGmshGrid(
    "triangular", "../../examples/meshes/sphere-h-0.1.msh")

# Create a space of piecewise constant basis functions over the grid.

pwiseConstants = lib.createPiecewiseConstantScalarSpace(context, grid)
pwiseLinears = lib.createPiecewiseLinearContinuousScalarSpace(context, grid)
pwiseLinears = pwiseConstants
#pwiseConstants = pwiseLinears

dirichletData = lib.createGridFunction(
    context, pwiseLinears, pwiseLinears, evalDirichletData)

# We now initialize the boundary operators.
# A boundary operator always takes three space arguments: a domain space,
# a range space, and the test space (dual to the range).
# Here, we just use L^2 projections. Hence, all spaces are identical.

slpOp = lib.createHelmholtz3dSingleLayerBoundaryOperator(
   context, pwiseConstants, pwiseLinears, pwiseConstants, k)
dlpOp = lib.createHelmholtz3dDoubleLayerBoundaryOperator(
    context, pwiseLinears, pwiseLinears, pwiseConstants, k)
adlpOp = lib.createHelmholtz3dAdjointDoubleLayerBoundaryOperator(
    context, pwiseLinears, pwiseLinears, pwiseConstants, k)
hypOp = lib.createHelmholtz3dHypersingularBoundaryOperator(
    context, pwiseLinears, pwiseLinears, pwiseConstants, k)
idOp = lib.createIdentityOperator(
    context, pwiseLinears, pwiseLinears, pwiseConstants)

# Standard arithmetic operators can be used to create linear combinations of
# boundary operators.

#lhsOp = idOp + 2 * adlpOp - 2j * k * slpOp
lhsOp = dlpOp

res = lhsOp*dirichletData

if 0:
	options = lib.createAssemblyOptions()
	options.switchToDenseMode()

	# The context object combines these settings

	context = lib.createContext(quadStrategy, options)

	slpOp = lib.createHelmholtz3dSingleLayerBoundaryOperator(
	   context, pwiseConstants, pwiseLinears, pwiseConstants, k)
	dlpOp = lib.createHelmholtz3dDoubleLayerBoundaryOperator(
	    context, pwiseLinears, pwiseLinears, pwiseConstants, k)
	adlpOp = lib.createHelmholtz3dAdjointDoubleLayerBoundaryOperator(
	    context, pwiseLinears, pwiseLinears, pwiseConstants, k)
	idOp = lib.createIdentityOperator(
	    context, pwiseLinears, pwiseLinears, pwiseConstants)
	lhsOp = dlpOp

	resaca = lhsOp*dirichletData
	res -= resaca
# Plot data

from bempp import visualization2 as vis
sphere = vis.tvtkGridFunction(res)
vis.plotScalarData(tvtkGridFunctions=sphere)
sys.exit()

#print dir(res.coefficients)
import matplotlib.pyplot as plt
plt.plot(abs(res.coefficients()))
plt.show()

sys.exit()




# Use the rhsData() Python function defined earlier to initialize the grid
# function that represents the right-hand side. The spaces are the domain space
# and the test space (in this case they are identical). rhsData() takes the
# surface normal as a parameter, so we set surfaceNormalDependent to True.

fun = lib.createGridFunction(
    context, pwiseConstants, pwiseConstants, rhsData, surfaceNormalDependent=True)


# We will now use GMRES to solve the problem.

# The default iterative solver supports several Krylov space methods.

solver = lib.createDefaultIterativeSolver(lhsOp)

# Create an initialization list for GMRES with tolerance 1e-5.
# A CG parameter list is also available for symmetric problems.

params = lib.defaultGmresParameterList(1e-5)
solver.initializeSolver(params)

# Solve...

solution = solver.solve(fun)
print solution.solverMessage() # Trilinos solver summary message

# ... and extract the solution. It is the normal derivative of the total field.

solfun = solution.gridFunction()


# We now do a simple 3d visualization. The visualization module provides some
# helper routines using the TVTK package. A lot more is possible by directly
# accessing the Mayavi visualization pipelines.
# (This example still uses the "old" visualization module. It needs to be
# updated to the "new" visualization module visualization2.)

# In order to evaluate the field we need to create a single-layer potential
# operator. Note: A single-layer potential operator is different from a
# single-layer boundary operator in the sense that its purpose is evaluation in
# free space, while the boundary operator is a weak form that lives on the
# boundary.

potential = lib.createHelmholtz3dSingleLayerPotentialOperator(context, k)

# We will create a three-plane view. For this we need to define the extents of
# the planes and the number of grid points in each dimension.

limits = (-5, 5)
dims = 200

# The total field is evaluated as u = u_{inc} - S u_{n}, where u_{inc} is the
# incident field and u_{n} the computed normal derivative stored in solfun.
# Hence, to get the total field we need to take the output of val = S u_{n} and
# compute u_{inc} - val for each grid point. This is done by the following
# 'transformation' function, which is given to the visualization routine.

def transformation(point, val):
    return np.real(np.exp(1j *k * point[:,0]) - val)

# We now plot the solution. The colorRange limits the color scale.

vis.plotThreePlanes(potential, solfun, limits, dims,
                    transformation=transformation, colorRange=(-1, 1))


