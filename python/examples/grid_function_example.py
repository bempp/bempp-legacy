import sys
sys.path.append("..")
import bempp
import math

# Define some functor

def fun(point):
    return math.cos(4 * point[0]) * math.sin(5 * point[1])

print "Importing grid..."
grid = bempp.GridFactory.importGmshGrid(
    "triangular", "../../examples/meshes/sphere-41440.msh")
space = bempp.piecewiseConstantScalarSpace(grid)
space.assignDofs()
ops = bempp.AssemblyOptions()

factory = bempp.defaultLocalAssemblerFactoryForOperatorsOnSurfaces()
f = bempp.gridFunctionFromSurfaceNormalIndependentFunction(
    space, space, fun, factory, ops)

f.exportToVtk(bempp.VtkWriter.CELL_DATA, "py_grid_fun", "py_gridfun")

