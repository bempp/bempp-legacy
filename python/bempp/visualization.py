# Copyright (C) 2011-2012 by the BEM++ Authors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import numpy as np

try:
    from tvtk.api import tvtk
    from mayavi import mlab
except ImportError:
    print "You need to have Enthought tvtk and mayavi installed for this module to work!"

def getTvtkGrid(grid):
    """Return a TVTK Object containing the grid"""

    if grid.topology=="triangular":
        (points,elems,auxData) = grid.leafView().getRawElementData()
        elem_list = elems[:-1,:].T
        mesh = tvtk.PolyData()
        mesh.points = points.T
        mesh.polys = elem_list
    else:
        raise TypeError("Visualization of this grid topology not implemented!")
    return mesh

def plotGrid(grid):
    """Plot a grid using Tvtk"""

    mesh = getTvtkGrid(grid)
    mapper = tvtk.DataSetMapper(input=mesh)
    actor  = tvtk.Actor(mapper=mapper)
    actor.property.representation = 'w'
    actor.property.ambient = 1
    v = mlab.figure()
    v.scene.add_actor(actor)
    mlab.show()
    
    
    
    

    
        
        

    
