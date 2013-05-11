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
from bempp import tools

try:
    from tvtk.api import tvtk
    from mayavi import mlab

    from traits.api import HasTraits, on_trait_change, Enum, Instance, List, Str, Bool, CBool, CFloat
    from traitsui.api import View, Item, SetEditor, Group
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


def plotGrid(grid,representation='wireframe'):
    """
    Plot a grid.
    """

    gridTvtkData = getTvtkGrid(grid)
    mlab.pipeline.surface(gridTvtkData,representation=representation)


def plotThreePlanes(points,vals):
    """plotThreePlanes(points,vals)

       Plot a three planes view of data.

       Parameters:
       -----------
       points       : An mgrid object defining points in a box
       vals         : The corresponding values

       Potential values in the correct format are returned by the function
       tools.evaluatePotentialInBox

    """



    from mayavi.tools.pipeline import scalar_field
    s = scalar_field(points[0],points[1],points[2],vals)

    mlab.pipeline.image_plane_widget(s,
                            plane_orientation='x_axes')
    mlab.pipeline.image_plane_widget(s,
                            plane_orientation='y_axes')
    mlab.pipeline.image_plane_widget(s,
                            plane_orientation='z_axes')
    mlab.outline()

