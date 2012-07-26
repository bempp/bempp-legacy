# Test script to evaluate a Helmholtz potential

import numpy as np
from tvtk.api import tvtk
from bempp import lib
from bempp import py_extensions as py_ext
from bempp import visualization as vis

k=2

def gridfundata(point,normal):
    return 2j*k*np.exp(k*1j*point[0])*(normal[0]-1)

# Create a potential

accuracy_options = lib.createAccuracyOptions()
accuracy_options.doubleRegular.orderIncrement=2
accuracy_options.doubleSingular.orderIncrement=1

factory = lib.createNumericalQuadratureStrategy("float64", "complex128",accuracy_options)
options = lib.createAssemblyOptions()
options.switchToAca(lib.core.AcaOptions())
context = lib.createContext(factory,options)
potential = lib.createHelmholtz3dSingleLayerPotentialOperator(context,k)


# Create the space

grid_factory = lib.createGridFactory()
grid = grid_factory.importGmshGrid("triangular","/Users/betcke/development/bempp/examples/meshes/sphere-614.msh")


pwiseConstants = lib.createPiecewiseConstantScalarSpace(context, grid)
pwiseConstants.assignDofs()
pwiseLinears = lib.createPiecewiseLinearContinuousScalarSpace(context, grid)
pwiseLinears.assignDofs()


# Create the single layer potential operator

slpOp = lib.createHelmholtz3dSingleLayerBoundaryOperator(context, pwiseConstants, pwiseConstants, pwiseConstants,k)
adlpOp = lib.createHelmholtz3dAdjointDoubleLayerBoundaryOperator(context, pwiseConstants, pwiseConstants, pwiseConstants, k)
id    = lib.createIdentityOperator(context, pwiseConstants, pwiseConstants, pwiseConstants)

# Create a grid function

fun = lib.gridFunctionFromSurfaceNormalDependentFunction(context, pwiseConstants, pwiseConstants, gridfundata)

# Create the operators

lhsOp = id+2*adlpOp-2j*k*slpOp


# Solve the problem

solver = lib.createDefaultIterativeSolver(lhsOp)
params = lib.defaultGmresParameterList(1e-5)
solver.initializeSolver(params)
solution = solver.solve(fun)
solfun = solution.gridFunction()

 

# Visualize

limits1 = (-5,5,-5,5,0,0)
limits2 = (0,0,-5,5,-5,5)
dims = (200,200)

(points1,vals1) = py_ext.evaluatePotentialOnPlane(potential,solfun,limits1,dims)
(points2,vals2) = py_ext.evaluatePotentialOnPlane(potential,solfun,limits2,dims)

vals1 = np.exp(1j*k*points1[:,0])-vals1
vals2 = np.exp(1j*k*points2[:,0])-vals2

# Now create a structured points object

g1 = tvtk.StructuredGrid(dimensions=(dims[0],dims[1],1),points=points1)
g2 = tvtk.StructuredGrid(dimensions=(dims[0],dims[1],1),points=points2)

# Add some data

g1.point_data.scalars = np.real(vals1)
g2.point_data.scalars = np.real(vals2)

# Render everything

mapper1 = tvtk.DataSetMapper(input=g1)
mapper1.scalar_range = (-1,1) #g.point_data.scalars.range
actor1 = tvtk.Actor(mapper=mapper1)

mapper2 = tvtk.DataSetMapper(input=g2)
mapper2.scalar_range = (-1,1) #g.point_data.scalars.range
actor2 = tvtk.Actor(mapper=mapper2)


grid_actor = vis.gridActor(grid)
legend = vis.legendActor(actor1)

vis.plotTvtkActors([actor1,actor2,grid_actor,legend])
 
