import sys
sys.path.append("..")
import bempp
import numpy

# Define some functor

def fun(point):
    return numpy.cos(4*point[0])*numpy.sin(5*point[1])

print "Importing grid..."
grid = bempp.GridFactory.importGmshGrid("triangular", "../../examples/meshes/sphere-41440.msh")
space=bempp.piecewiseConstantScalarSpace(grid)
space.assignDofs()
ops=bempp.AssemblyOptions()

factory=bempp.standardLocalAssemblerFactoryForOperatorsOnSurfaces()
f=bempp.surfaceNormalIndependentFunctor(fun,'float64',3,1)
g=bempp.gridFunctionFromSurfaceNormalIndependentFunctor(space,f,factory,ops)

g.exportToVtk(bempp.VtkWriter.CELL_DATA,"py_grid_fun","py_gridfun")

