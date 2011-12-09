#!/usr/bin/env python

# This script imports a grid in the Gmsh format and exports the z coordinate 
# of each vertex to the file "output.vtu" in the VTK format.

# Help Python find the bempp module
import sys
sys.path.append("../python/")

import bempp
import numpy as np

# Create the grid, the leaf view and its index set
print "Importing grid..."
grid = bempp.GridFactory.importGmshGrid("triangular", "head.gmsh")
view = grid.leafView()
index_set = view.indexSet()

# Initialise the container for the data to be exported
n_vertices = view.entityCount(2)
data = np.zeros((1, n_vertices))

# To each vertex assign its y coordinate
print "Traversing grid..."
for e in view.entities(2):
    geo = e.geometry()
    ctr = geo.center()
    data[0, index_set.entityIndex(e)] = ctr[1] # the z coordinate

# Export data
vtk_writer = view.vtkWriter()
vtk_writer.addVertexData(data, "y_coord")
print "Exporting data..."
outName = vtk_writer.write("output")
print "Data exported to file '%s'." % outName

