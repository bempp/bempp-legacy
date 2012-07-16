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

    if grid.topology()=="triangular":
        (points,elems,auxData) = grid.leafView().getRawElementData()
        elem_list = elems[:-1,:].T
        mesh = tvtk.PolyData()
        mesh.points = points.T
        mesh.polys = elem_list
    else:
        raise TypeError("Visualization of this grid topology not implemented!")
    return mesh

def plotTvtkActors(tvtkActors):
    """Plot a number of TVTK actors in the same plot"""

    import collections

    v = mlab.figure()
    if isinstance(tvtkActors, collections.Iterable):
        for actor in tvtkActors: v.scene.add_actor(actor) # Input is iterable
    else:
        v.scene.add_actor(tvtkActors)  # Input is not iteratble
    mlab.show()


def gridActor(grid):
    """Return a grid actor using TVTK"""

    mesh = getTvtkGrid(grid)
    mapper = tvtk.DataSetMapper(input=mesh)
    actor  = tvtk.Actor(mapper=mapper)
    actor.property.representation = 'w'
    actor.property.ambient = 1
    return actor
    

def gridFunctionActor(gridFun,data_type='vertex_data',transformation='real'):
    """Plot a grid function usint TVTK"""

    if not data_type in ["cell_data", "vertex_data"]:
        raise ValueError("Unknown mode specified. Valid modes are 'vertex_data' and 'cell_data'!")

    if not hasattr(transformation, '__call__'):
        if transformation=='real':
            data_transform = lambda x:np.real(x)
        elif transformation=='imag':
            data_transform = lambda x:np.imag(x)
        elif transformation=='abs':
            data_transform = lambda x:np.abs(x)
        else:
            raise ValueError("Unknown value for 'transformation'. It needs to be 'real', 'imag', 'abs' or a Python Callable!")
    else:
        data_transform = transformation

    data_transform = np.vectorize(data_transform)

    mesh = getTvtkGrid(gridFun.grid())
    if data_type=="vertex_data":
        values = gridFun.evaluateAtSpecialPoints("vertex_data").flatten()
        tvtk_data = mesh.point_data
    elif data_type=="cell_data":
        values = gridFun.evaluateAtSpecialPoints("cell_data").flatten()
        tvtk_data = mesh.cell_data

    values = data_transform(values)

    tvtk_data.scalars = values
    mapper = tvtk.DataSetMapper(input = mesh)
    mapper.scalar_range = tvtk_data.scalars.range
    actor = tvtk.Actor(mapper=mapper)
    return actor

def legendActor(actor):
    """Return a legend object for the specified actor"""

    scalar_bar = tvtk.ScalarBarActor()
    scalar_bar.lookup_table = actor.mapper.lookup_table
    return scalar_bar

def plotGridFunction(*args,**kwargs):
    """Simple grid function plotter"""

    fun = gridFunctionActor(*args,**kwargs)
    legend = legendActor(fun)
    plotTvtkActors([fun,legend])
    
    

    

    
    
    
    

    
        
        

    
