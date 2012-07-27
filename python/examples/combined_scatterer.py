# Scattering from a sphere using a combined direct formulation

# Import numpy and library modules

import numpy as np 
from bempp import lib
from bempp import visualization as vis

# Define the wavenumber

k=6 

# The rhs of the combined formulation

def gridfundata(point,normal):
    return 2j*k*np.exp(k*1j*point[0])*(normal[0]-1)


# Set accuracy options

accuracy_options = lib.createAccuracyOptions()
accuracy_options.doubleRegular.orderIncrement=2 # 2 orders higher than default accuracy for regular integrals
accuracy_options.doubleSingular.orderIncrement=1 # 1 order higher than default accuracy for singular integrals


# The default factory for numerical quadrature
# "float64" is the basis function type and "complex128" is the result type"
 
factory = lib.createNumericalQuadratureStrategy("float64", "complex128",accuracy_options)

# We want ACA assembly

options = lib.createAssemblyOptions()
options.switchToAca(lib.createAcaOptions())

# The context object combines all information

context = lib.createContext(factory,options)


# Load a grid

grid_factory = lib.createGridFactory()
grid = grid_factory.importGmshGrid("triangular","../../examples/meshes/sphere-2590.msh")

# Create a space of piecewise constant basis functions over the grid
# Currently assignDofs() always needs to be invoked manually to assign the degrees of freedom in the space

pwiseConstants = lib.createPiecewiseConstantScalarSpace(context, grid)
pwiseConstants.assignDofs()

# We now initialize the potential operators and the identity operator.
# A potential always takes three space arguments, a domain space, a range space, and the test space (dual to the range).
# Here, we just use L^2 projections. Hence, all spaces are identical.

slpOp = lib.createHelmholtz3dSingleLayerBoundaryOperator(context, pwiseConstants, pwiseConstants, pwiseConstants,k)
adlpOp = lib.createHelmholtz3dAdjointDoubleLayerBoundaryOperator(context, pwiseConstants, pwiseConstants, pwiseConstants, k)
id    = lib.createIdentityOperator(context, pwiseConstants, pwiseConstants, pwiseConstants)

# The operators support simple arithmetic expressions to create linear combinations of operators.

lhsOp = id+2*adlpOp-2j*k*slpOp

# We now initialize the grid function that represents the right-hand side.
# In the combined formulation the rhs also needs the normal. Hence, we take the function that
# turns a Python callable into a SurfaceNormalDependent grid function.
# The spaces are the domain space and the test space (in this case they are identical).

fun = lib.gridFunctionFromSurfaceNormalDependentFunction(context, pwiseConstants, pwiseConstants, gridfundata)


# We now use GMRES to solve the problem

# The default iterative solver supports several Krylov space methods.

solver = lib.createDefaultIterativeSolver(lhsOp)

# Create an initialization list for Gmres with default tolerance 1E-5.
# A CG parameter list is also available for symmetric problems.

params = lib.defaultGmresParameterList(1e-5)
solver.initializeSolver(params)

# Now solve and extract the solution

solution = solver.solve(fun)
print solution.solverMessage() # Trilinos solver summary message

# Now extract the solution. It is the normal derivative of the total field.
 
solfun = solution.gridFunction()


# We now do a simple 3d visualization.
# Some routines are implemented to access the TVTK visualization routines.
# A lot more is possible by directly accessing the Mayavi visualization pipelines.

# In order to evaluate the field we need to create a single layer potential.
# Note: A single layer potential is different from a single layer boundary operator
# in the sense that its purpose is evaluation in free space, while the boundary
# operator is a weak form that lives on the boundary.

potential = lib.createHelmholtz3dSingleLayerPotentialOperator(context,k)

# We will create a three plane view. For this we need to define the extents of the planes
# and the number of grid points in each dimension.

limits = (-5,5,-5,5)
dims = (200,200)

# The total field is evaluated as u=u_{inc}-Su_{n}, where u_{inc} is the incident
# field and u_{n} the computed normal derivative stored in solfun. Hence, to
# get the total field we need to take the output of val=Su_{n} and compute
# u_{inc}-val for each grid point. This is done by the following 'transformation'
# function, which is given to the visualization routine.

def transformation(point,val):
    return np.real(np.exp(1j*k*point[:,0])-val)

# We now plot the solution. The colorRange limits the color scale. 

vis.plotThreePlanes(potential,solfun,limits,dims,transformation=transformation,
                    colorRange=(-1,1))

 
