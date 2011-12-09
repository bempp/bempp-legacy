#!/usr/bin/env python

# This script creates a structured triangular grid covering a unit square
# and exports the y coordinate of each vertex in to the file "output.vtu" 
# in the VTK format.

# Help Python find the bempp module
import sys
sys.path.append("../python/")

import bempp
import numpy as np

# Create the grid, the leaf view and its index set
grid = bempp.GridFactory.createStructuredGrid("triangular", (0, 0), (1, 1), (5, 6))
view = grid.leafView()
index_set = view.indexSet()

# Initialise the container for the data to be exported
n_vertices = view.entityCount(2)
data = np.zeros((1, n_vertices))

# To each vertex assign its y coordinate
for e in view.entities(2):
    geo = e.geometry()
    ctr = geo.center()
    data[0, index_set.entityIndex(e)] = ctr[1] # the y coordinate
print "Data:\n", data

# Export data
vtk_writer = view.vtkWriter()
vtk_writer.addVertexData(data, "y_coord")
outName = vtk_writer.write("output")
print "Data exported to file '%s'." % outName

