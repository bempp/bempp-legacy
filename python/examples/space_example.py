#!/usr/bin/env python

# This script imports a grid in the Gmsh format and exports the z coordinate 
# of each vertex to the file "output.vtu" in the VTK format.

# Help Python find the bempp module
import sys
sys.path.append("..")

import bempp
import numpy as np

# Create the grid, the leaf view and its index set
print "Importing grid..."
grid = bempp.GridFactory.importGmshGrid("triangular", "../../examples/meshes/sphere-614.msh")
space=bempp.PiecewiseConstantScalarSpace('float32',grid)
space.assignDofs()
print space.dofsAssigned()
print space.templateType
print object.__getattribute__(space,'impl')


