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
import bempp.py_extensions as py_ext

try:
    from tvtk.api import tvtk
    from mayavi import mlab
except ImportError:
    print "You need to have Enthought tvtk and mayavi installed for this module to work!"

def getTvtkGrid(grid):
    """Return a TVTK object containing the grid"""

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
    """Plot a number of TVTK actors in the same plot."""

    import collections

    v = mlab.figure()
    if isinstance(tvtkActors, collections.Iterable):
        for actor in tvtkActors: v.scene.add_actor(actor) # Input is iterable
    else:
        v.scene.add_actor(tvtkActors)  # Input is not iteratble
    mlab.show()


def gridActor(grid):
    """
    Return a TVTK actor representing a grid.

    *Parameters:*
       - grid (Grid)
            A BEM++ grid object.
    """

    mesh = getTvtkGrid(grid)
    mapper = tvtk.DataSetMapper(input=mesh)
    actor  = tvtk.Actor(mapper=mapper)
    actor.property.representation = 'w'
    actor.property.ambient = 1
    return actor

def gridFunctionActor(gridFun,dataType='vertex_data',transformation='real'):
    """
    Return a TVTK actor representing a grid function.

    *Parameters:*
       - gridFun (gridFunction).
            The grid function to be plotted.
       - dataType ('vertex_data' or 'cell_data')
            Determines whether the plot should be constructed from the function
            values at element vertices or at element centres.
       - transformation ('real', 'imag', 'abs' or a callable object)
            Determines how the function is transformed before plotting.
    """

    if not dataType in ["cell_data", "vertex_data"]:
        raise ValueError("Unknown mode specified. Valid modes are 'vertex_data' and 'cell_data'!")

    if not hasattr(transformation, '__call__'):
        if transformation=='real':
            data_transform = lambda x:np.real(x)
        elif transformation=='imag':
            data_transform = lambda x:np.imag(x)
        elif transformation=='abs':
            data_transform = lambda x:np.abs(x)
        else:
            raise ValueError("Unknown value for 'transformation'. It needs to be 'real', 'imag', 'abs' or a Python callable!")
    else:
        data_transform = transformation

    mesh = getTvtkGrid(gridFun.grid())
    if dataType=="vertex_data":
        values = gridFun.evaluateAtSpecialPoints("vertex_data").flatten()
        tvtk_data = mesh.point_data
    elif dataType=="cell_data":
        values = gridFun.evaluateAtSpecialPoints("cell_data").flatten()
        tvtk_data = mesh.cell_data

    values = data_transform(values)

    tvtk_data.scalars = values
    mapper = tvtk.DataSetMapper(input = mesh)
    mapper.scalar_range = tvtk_data.scalars.range
    actor = tvtk.Actor(mapper=mapper)
    return actor

def scalarDataOnRegularGridActor(
        points, data, dimensions,
        colorRange=None,transformation='real'):
    """
    Return a TVTK actor representing the plot of a function interpolated on
    a regular grid.
    """

    if points.shape[0] != 3 or points.ndim != 2:
        raise ValueError("incorrect shape")
    data = data.squeeze()
    if data.ndim != 1:
        raise ValueError("incorrect shape")

    if not hasattr(transformation, '__call__'):
        if transformation=='real':
            data_transform = lambda point,val:np.real(val)
        elif transformation=='imag':
            data_transform = lambda point,val:np.imag(val)
        elif transformation=='abs':
            data_transform = lambda point,val:np.abs(val)
        else:
            raise ValueError("Unknown value for 'transformation'. It needs to be 'real', 'imag', 'abs' or a Python Callable!")
    else:
        data_transform = transformation

    data = data_transform(points,data)
    dims = dimensions

    if colorRange is None:
        minVal = np.min(data)
        maxVal = np.max(data)
        colorRange = (minVal, maxVal)

    g = tvtk.StructuredGrid(dimensions=(dims[1], dims[0], 1), points=points.T)
    g.point_data.scalars = data

    # Create actors
    mapper = tvtk.DataSetMapper(input=g)
    mapper.scalar_range = colorRange
    return tvtk.Actor(mapper=mapper)

def legendActor(actor):
    """Return a TVTK actor representing the legend of another actor."""

    scalar_bar = tvtk.ScalarBarActor()
    scalar_bar.lookup_table = actor.mapper.lookup_table
    return scalar_bar

def plotgridFunction(*args,**kwargs):
    """
    Plot a grid function.

    This function takes the same parameters as gridFunctionActor().
    """

    fun = gridFunctionActor(*args,**kwargs)
    legend = legendActor(fun)
    plotTvtkActors([fun,legend])
    
def plotGrid(grid):
    """
    Plot a grid.

    This function takes the same parameters as gridActor().
    """

    grid_actor = gridActor(grid)
    plotTvtkActors(grid_actor)

def plotThreePlanes(potentialOp, gridFun, limits, dimensions,
                    colorRange=None, transformation='real', evalOps=None):
    """
    Plot the potential generated by applying a potential operator to a grid
    function on the xy, xz and yz planes.

    *Parameters:*
       - potentialOp (PotentialOperator)
            A potential operator.
       - gridFun (GridFunction)
            A grid function.
       - limits (tuple)
            Tuple (min, max) or (xmin, xmax, ymin, ymax, zmin, zmax)
            specifying the extent of each plane on which the potential
            will be plotted.
       - dimensions (tuple)
            Scalar or tuple (xdim, ydim, zdim) specifying the number of samples
            in each direction.
       - colorRange (tuple)
            Tuple (min, max) determining the range of data to be plotted.
            If set to None, the data range is determined automatically.
       - transformation ('real', 'imag', 'abs' or a callable object)
            Determines how the potential is transformed before plotting.
       - evalOps (EvaluationOptions)
            Options controlling the evaluation of the potential.
    """

    if np.isscalar(dimensions):
        dims = (dimensions, dimensions, dimensions)
    else:
        if len(dimensions) == 3:
            dims = dimensions
        else:
            raise ValueError("dimensions must be a scalar or a tuple with 3 elements")
    if len(limits) == 2:
        lims = (limits[0], limits[1], limits[0], limits[1], limits[0], limits[1])
    elif len(limits) == 6:
        lims = limits
    else:
        raise ValueError("limits must be a tuple with 2 or 6 elements")
    origin = ((lims[0] + lims[1]) / 2.,
              (lims[2] + lims[3]) / 2.,
              (lims[4] + lims[5]) / 2.)
    (points1,vals1) = py_ext.evaluatePotentialOnPlane(
        potentialOp,gridFun,lims[:4],dims[:2],plane="xy",origin=origin,
        evalOps=evalOps)
    (points2,vals2) = py_ext.evaluatePotentialOnPlane(
        potentialOp,gridFun,lims[:2]+lims[4:],dims[:1]+dims[2:],plane="xz",
        origin=origin,evalOps=evalOps)
    (points3,vals3) = py_ext.evaluatePotentialOnPlane(
        potentialOp,gridFun,lims[2:],dims[1:],plane="yz",origin=origin,
        evalOps=evalOps)

    if not hasattr(transformation, '__call__'):
        if transformation=='real':
            data_transform = lambda point,val:np.real(val)
        elif transformation=='imag':
            data_transform = lambda point,val:np.imag(val)
        elif transformation=='abs':
            data_transform = lambda point,val:np.abs(val)
        else:
            raise ValueError("Unknown value for 'transformation'. It needs to be 'real', 'imag', 'abs' or a Python Callable!")
    else:
        data_transform = transformation

    vals1 = data_transform(points1,vals1)
    vals2 = data_transform(points2,vals2)
    vals3 = data_transform(points3,vals3)

    if colorRange is None:
        minVal = np.min([vals1,vals2,vals3])
        maxVal = np.max([vals1,vals2,vals3])
        colorRange = (minVal,maxVal)

    g1 = tvtk.StructuredGrid(dimensions=(dims[0],dims[1],1),points=points1)
    g2 = tvtk.StructuredGrid(dimensions=(dims[0],dims[2],1),points=points2)
    g3 = tvtk.StructuredGrid(dimensions=(dims[1],dims[2],1),points=points3)


    # Add data

    g1.point_data.scalars = vals1
    g2.point_data.scalars = vals2
    g3.point_data.scalars = vals3


    # Create actors

    mapper1 = tvtk.DataSetMapper(input=g1)
    mapper1.scalar_range = colorRange
    actor1 = tvtk.Actor(mapper=mapper1)

    mapper2 = tvtk.DataSetMapper(input=g2)
    mapper2.scalar_range = colorRange
    actor2 = tvtk.Actor(mapper=mapper2)

    mapper3 = tvtk.DataSetMapper(input=g3)
    mapper3.scalar_range = colorRange 
    actor3 = tvtk.Actor(mapper=mapper3)

    gActor = gridActor(gridFun.grid())

    legend = legendActor(actor1)

    plotTvtkActors([actor1,actor2,actor3,gActor,legend])
