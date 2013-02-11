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

    from traits.api import HasTraits, on_trait_change, Enum, Instance, List, Str, Bool, CBool, CFloat
    from traitsui.api import View, Item, SetEditor, Group
    from traitsui.wx.check_list_editor import SimpleEditor
    from tvtk.pyface.scene_editor import SceneEditor
    from mayavi.tools.mlab_scene_model import MlabSceneModel
    from mayavi.core.ui.mayavi_scene import MayaviScene
    from mayavi.sources.api import VTKDataSource

except ImportError:
    print "You need to have Enthought tvtk, mayavi and traits installed for this module to work!"

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

class _VectorVisualization(HasTraits):

    real_imag = Enum('Real Part of Vector Field','Imaginary Part of Vector Field')
    legend = Enum('Scalar Legend','Vector Legend')
    point_cell = Enum('Point Data','Cell Data')
    scene = Instance(MlabSceneModel, ())
    enable_surface = Bool(True)
    enable_vectors = Bool(False)
    enable_scalars = Bool(True)
    enable_grid    = Bool(False)
    vector_scale_size = CFloat(0.1)

    def __init__(self,g):
        HasTraits.__init__(self)
        self.src = VTKDataSource(data=g)

    @on_trait_change('scene.activated')
    def create_plot(self):
        from mayavi.modules.api import Surface, Vectors
        self.surface = Surface()
        self.surface1 = Surface()
        self.surface1.actor.property.representation = 'wireframe'
        self.surface1.actor.actor.visibility = self.enable_grid
        self.surface.actor.actor.visibility=self.enable_surface
        self.vectors = Vectors()
        self.engine = self.scene.engine
        self.engine.add_source(self.src)
        self.engine.add_module(self.surface, obj=self.src)
        self.engine.add_module(self.surface1, obj=self.src)
        self.engine.add_module(self.vectors, obj=self.src)
        self.module_manager = self.engine.scenes[0].children[0].children[0]
        self.module_manager.vector_lut_manager.show_legend = True
        self.vectors.actor.actor.visibility=self.enable_vectors
        self.surface.actor.actor.visibility=self.enable_surface
        if self.legend == "Scalar Legend":
            self.module_manager.vector_lut_manager.show_legend = False
            self.module_manager.scalar_lut_manager.show_legend = True
        else:
            self.module_manager.vector_lut_manager.show_legend = True
            self.module_manager.scalar_lut_manager.show_legend = False
        if self.point_cell == "Point Data":
            self.module_manager.lut_data_mode = 'point data'
        else:
            self.module_manager.lut_data_mode = 'cell data'
        self.vectors.glyph.glyph.scale_factor = self.vector_scale_size
        self.src.point_scalars_name = 'abs^2'
        self.src.cell_scalars_name = 'abs^2'


    @on_trait_change('real_imag')
    def update_real_imag(self):
        if self.real_imag=="Real Part of Vector Field":
            self.src.point_vectors_name = 'real'
            self.src.cell_vectors_name = 'real'
        else:
            self.src.point_vectors_name = 'imag'
            self.src.cell_vectors_name = 'imag'


    @on_trait_change('enable_grid')
    def update_grid(self):
        self.surface1.actor.actor.visibility = self.enable_grid

    @on_trait_change('point_cell')
    def update_point_cell(self):
        if self.point_cell == "Point Data":
            self.module_manager.lut_data_mode = 'point data'
        else:
            self.module_manager.lut_data_mode = 'cell data'


    @on_trait_change('enable_surface')
    def update_surface(self):
        self.surface.actor.actor.visibility=self.enable_surface

    @on_trait_change('enable_vectors')
    def update_vectors(self):
        self.vectors.actor.actor.visibility=self.enable_vectors
        if self.enable_vectors:
            self.legend = "Vector Legend"
        else:
            self.legend = "Scalar Legend"

    @on_trait_change('legend')
    def update_legend(self):
        if self.legend == "Scalar Legend":
            self.module_manager.vector_lut_manager.show_legend = False
            self.module_manager.scalar_lut_manager.show_legend = True
        else:
            self.enable_vectors = True
            self.module_manager.vector_lut_manager.show_legend = True
            self.module_manager.scalar_lut_manager.show_legend = False

    @on_trait_change('vector_scale_size')
    def update_vector_scale_size(self):
        if self.vector_scale_size>0:
            self.vectors.glyph.glyph.scale_factor = self.vector_scale_size


    view = View(Item(name='scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=500, width=500, show_label=False),
                Group(
                    Item(name="real_imag",style='custom',show_label=False),
                    Item(name="legend",style='custom',show_label=False),
                Group(Item(name="point_cell",show_label=False),
                      Item(name="vector_scale_size",label="Vector Scale"),
                Group(
                    Item(name="enable_surface",label="Display scalar density"),
                    Item(name="enable_vectors",label="Enable vectors"),
                    Item(name='enable_grid',label="Show grid"),orientation="horizontal")
),
                orientation="vertical"),
                resizable=True,title="Grid Function Viewer")


class _ScalarVisualization(HasTraits):

    real_imag = Enum('Real Part','Imaginary Part', 'Square Density')
    enable_legend = Bool(True)
    point_cell = Enum('Point Data','Cell Data')
    scene = Instance(MlabSceneModel, ())
    enable_surface = Bool(True)
    enable_scalars = Bool(True)
    enable_grid    = Bool(False)

    def __init__(self,g):
        HasTraits.__init__(self)
        self.src = VTKDataSource(data=g)

    @on_trait_change('scene.activated')
    def create_plot(self):
        from mayavi.modules.api import Surface, Vectors
        self.surface = Surface()
        self.surface1 = Surface()
        self.surface1.actor.property.representation = 'wireframe'
        self.surface1.actor.actor.visibility = self.enable_grid
        self.surface.actor.actor.visibility=self.enable_surface
        self.engine = self.scene.engine
        self.engine.add_source(self.src)
        self.engine.add_module(self.surface, obj=self.src)
        self.engine.add_module(self.surface1, obj=self.src)
        self.module_manager = self.engine.scenes[0].children[0].children[0]
        self.surface.actor.actor.visibility=self.enable_surface
        self.module_manager.scalar_lut_manager.show_legend = self.enable_legend
        if self.point_cell == "Point Data":
            self.module_manager.lut_data_mode = 'point data'
        else:
            self.module_manager.lut_data_mode = 'cell data'


    @on_trait_change('real_imag')
    def update_real_imag(self):
        if self.real_imag=="Real Part":
            self.src.point_scalars_name = 'real'
            self.src.cell_scalars_name = 'real'
        elif self.real_imag=="Imaginary Part":
            self.src.point_scalars_name = 'imag'
            self.src.cell_scalars_name = 'imag'
        else:
            self.src.point_scalars_name = 'abs^2'
            self.src.cell_scalars_name = 'abs^2'


    @on_trait_change('enable_grid')
    def update_grid(self):
        self.surface1.actor.actor.visibility = self.enable_grid

    @on_trait_change('point_cell')
    def update_point_cell(self):
        if self.point_cell == "Point Data":
            self.module_manager.lut_data_mode = 'point data'
        else:
            self.module_manager.lut_data_mode = 'cell data'


    @on_trait_change('enable_surface')
    def update_surface(self):
        self.surface.actor.actor.visibility=self.enable_surface

    @on_trait_change('enable_legend')
    def update_legend(self):
        self.module_manager.scalar_lut_manager.show_legend = self.enable_legend


    view = View(Item(name='scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=500, width=500, show_label=False),
                Group(
                    Item(name="real_imag",style='custom',show_label=False),
                Group(Item(name="point_cell",show_label=False),
                Group(
                    Item(name="enable_surface",label="Display scalar density"),
                    Item(name='enable_grid',label="Show grid"),
                    Item(name="enable_legend", label = "Show legend" ),
                    orientation="horizontal")
),
                orientation="vertical"),
                resizable=True,title="Grid Function Viewer")

def plotGridFunction(g):
    """Visualize a Grid Function. Returns a Traits object that contains the visualization."""


    tvtkObj = getTvtkGrid(g.grid())
    point_data = g.evaluateAtSpecialPoints("vertex_data")
    cell_data = g.evaluateAtSpecialPoints("cell_data")

    tvtkObj.cell_data.add_array(np.real(cell_data.T))
    tvtkObj.cell_data.add_array(np.imag(cell_data.T))
    tvtkObj.cell_data.add_array(np.sum(abs(cell_data)**2,axis=0))
    tvtkObj.cell_data.get_abstract_array(0).name = 'real'
    tvtkObj.cell_data.get_abstract_array(1).name = 'imag'
    tvtkObj.cell_data.get_abstract_array(2).name = 'abs^2'

    tvtkObj.point_data.add_array(np.real(point_data.T))
    tvtkObj.point_data.add_array(np.imag(point_data.T))
    tvtkObj.point_data.add_array(np.sum(abs(point_data)**2,axis=0))
    tvtkObj.point_data.get_abstract_array(0).name = 'real'
    tvtkObj.point_data.get_abstract_array(1).name = 'imag'
    tvtkObj.point_data.get_abstract_array(2).name = 'abs^2'

    if g.componentCount()==3:
        tvtkObj.cell_data.set_active_scalars('abs^2')
        tvtkObj.point_data.set_active_scalars('abs^2')
        tvtkObj.cell_data.set_active_vectors('real')
        tvtkObj.point_data.set_active_vectors('imag')
        v = _VectorVisualization(tvtkObj)

    elif g.componentCount()==1:
        tvtkObj.cell_data.set_active_scalars('real')
        tvtkObj.point_data.set_active_scalars('real')
        v =  _ScalarVisualization(tvtkObj)
    else:
        raise Exception("plotGridFunction: Only GridFunctions with componentCount 1 or 3 are supported.")

    v.configure_traits()
    return v


def plotGrid(grid):
    """
    Plot a grid.
    """

    gridTvtkData = getTvtkGrid(grid)
    mlab.pipeline.surface(gridTvtkData,representation='wireframe')

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
