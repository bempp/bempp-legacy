#!/usr/bin/env python

# This script imports a grid in the Gmsh format and exports the z coordinate 
# of each vertex to the file "output.vtu" in the VTK format.

# Help Python find the bempp module
import sys
sys.path.append("..")

import bempp
import numpy as np

print "Importing grid..."
grid = bempp.GridFactory.importGmshGrid("triangular", "../../examples/meshes/sphere-614.msh")
space = bempp.piecewiseConstantScalarSpace(grid)
space.assignDofs()
print space.dofsAssigned()
#print space._basisFunctionType

op = bempp.laplace3dSingleLayerPotentialOperator(space, space, space)
print "operator has a singular kernel?", not op.isRegular()

factory = bempp.standardLocalAssemblerFactoryForOperatorsOnSurfaces()
options = bempp.AssemblyOptions()
options.switchToAca(bempp.AcaOptions())
op.assembleWeakForm(factory, options)
wf = op.weakForm()
