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

    from traits.api import HasTraits, on_trait_change, Enum, Instance, List, Str, Bool, CBool, CFloat, TraitError
    from traitsui.api import View, Item, SetEditor, Group
    from tvtk.pyface.scene_editor import SceneEditor
    from mayavi.tools.mlab_scene_model import MlabSceneModel
    from mayavi.core.ui.mayavi_scene import MayaviScene
    from mayavi.sources.api import VTKDataSource
except ImportError:
    print "You need to have Enthought tvtk, mayavi and traits installed for this module to work!"

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

    def __init__(self, tvtkGridSrcs=None, tvtkGridFunctionSrcs=None,
                 tvtkStructuredGridDataSrcs=None, dataRange=None):
        HasTraits.__init__(self)
        if tvtkGridSrcs is not None:
            try:
                self.tvtkGridSrcs = [VTKDataSource(data=tvtkSrc) for tvtkSrc in tvtkGridSrcs]
            except TypeError:
                self.tvtkGridSrcs = [VTKDataSource(data=tvtkGridSrcs)]
        else:
            self.tvtkGridSrcs = []
        if tvtkGridFunctionSrcs is not None:
            try:
                self.tvtkGridFunctionSrcs = [VTKDataSource(data=tvtkSrc)
                                             for tvtkSrc in tvtkGridFunctionSrcs]
            except TypeError:
                self.tvtkGridFunctionSrcs = [VTKDataSource(data=tvtkGridFunctionSrcs)]
        else:
            self.tvtkGridFunctionSrcs = []
        if tvtkStructuredGridDataSrcs is not None:
            try:
                self.tvtkStructuredGridDataSrcs = [VTKDataSource(data=tvtkSrc)
                                             for tvtkSrc in tvtkStructuredGridDataSrcs]
            except TypeError:
                self.tvtkStructuredGridDataSrcs = [VTKDataSource(data=tvtkStructuredGridDataSrcs)]
        else:
            self.tvtkStructuredGridDataSrcs = []
        self.dataRange = dataRange

    @on_trait_change('scene.activated')
    def create_plot(self):
        from mayavi.modules.api import Surface, Vectors
        self.engine = self.scene.engine
        self.surface = []
        self.surface1 = []
        self.vectors = []
        # self.surface = Surface()
        # self.surface1 = Surface()
        # self.surface1.actor.property.representation = 'wireframe'
        # self.surface1.actor.actor.visibility = self.enable_grid
        # self.surface.actor.actor.visibility=self.enable_surface
        # self.vectors = Vectors()
        # self.engine = self.scene.engine
        # self.engine.add_source(self.src)
        # self.engine.add_module(self.surface, obj=self.src)
        # self.engine.add_module(self.surface1, obj=self.src)
        # self.engine.add_module(self.vectors, obj=self.src)
        # self.module_manager = self.engine.scenes[0].children[0].children[0]
        # self.module_manager.vector_lut_manager.show_legend = True
        # self.vectors.actor.actor.visibility=self.enable_vectors
        # self.surface.actor.actor.visibility=self.enable_surface
        for src in self.tvtkGridFunctionSrcs + self.tvtkStructuredGridDataSrcs:
            self.surface.append(Surface())
            self.vectors.append(Vectors())
            self.surface[-1].actor.actor.visibility=self.enable_surface
            self.vectors[-1].actor.actor.visibility=self.enable_vectors
            self.vectors[-1].glyph.glyph.scale_factor = self.vector_scale_size
            self.engine.add_source(src)
            self.engine.add_module(self.surface[-1], obj=src)
            self.engine.add_module(self.vectors[-1], obj=src)
        for src in self.tvtkGridSrcs:
            self.surface1.append(Surface())
            self.surface1[-1].actor.property.representation = 'wireframe'
            self.surface1[-1].actor.actor.visibility = self.enable_grid
            self.surface1[-1].actor.mapper.scalar_visibility = False
            self.engine.add_source(src)
            self.engine.add_module(self.surface1[-1], obj=src)
        self.module_manager = self.engine.scenes[0].children[0].children[0]
        if self.dataRange is not None:
            self.module_manager.vector_lut_manager.data_range = self.dataRange
            self.module_manager.scalar_lut_manager.data_range = self.dataRange
            self.module_manager.vector_lut_manager.use_default_range = False
            self.module_manager.scalar_lut_manager.use_default_range = False
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

    @on_trait_change('real_imag')
    def update_real_imag(self):
        if self.real_imag=="Real Part of Vector Field":
            for src in self.tvtkGridFunctionSrcs + self.tvtkStructuredGridDataSrcs:
                src.point_vectors_name = 'real'
                try:
                    src.cell_vectors_name = 'real'
                except TraitError:
                    pass # cell data not present
        else:
            for src in self.tvtkGridFunctionSrcs + self.tvtkStructuredGridDataSrcs:
                src.point_vectors_name = 'imag'
                try:
                    src.cell_vectors_name = 'imag'
                except TraitError:
                    pass # cell data not present

    @on_trait_change('enable_grid')
    def update_grid(self):
        for s in self.surface1:
            s.actor.actor.visibility = self.enable_grid

    @on_trait_change('point_cell')
    def update_point_cell(self):
        if self.point_cell == "Point Data":
            self.module_manager.lut_data_mode = 'point data'
        else:
            self.module_manager.lut_data_mode = 'cell data'

    @on_trait_change('enable_surface')
    def update_surface(self):
        for s in self.surface:
            s.actor.actor.visibility = self.enable_surface

    @on_trait_change('enable_vectors')
    def update_vectors(self):
        for s in self.vectors:
            s.actor.actor.visibility = self.enable_vectors
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
            for s in self.vectors:
                s.glyph.glyph.scale_factor = self.vector_scale_size

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
                    Item(name='enable_grid',label="Show grid"),orientation="horizontal")),
                orientation="vertical"),
                resizable=True,title="Grid Function Viewer")


class _ScalarVisualization(HasTraits):

    real_imag = Enum('Real Part', 'Imaginary Part', 'Squared Density')
    enable_legend = Bool(True)
    point_cell = Enum('Point Data', 'Cell Data')
    scene = Instance(MlabSceneModel, ())
    enable_surface = Bool(True)
    enable_scalars = Bool(True)
    enable_grid    = Bool(False)

    def __init__(self, tvtkGridSrcs=None, tvtkGridFunctionSrcs=None,
                 tvtkStructuredGridDataSrcs=None, dataRange=None):
        HasTraits.__init__(self)
        if tvtkGridSrcs is not None:
            try:
                self.tvtkGridSrcs = [VTKDataSource(data=tvtkSrc) for tvtkSrc in tvtkGridSrcs]
            except TypeError:
                self.tvtkGridSrcs = [VTKDataSource(data=tvtkGridSrcs)]
        else:
            self.tvtkGridSrcs = []
        if tvtkGridFunctionSrcs is not None:
            try:
                self.tvtkGridFunctionSrcs = [VTKDataSource(data=tvtkSrc)
                                             for tvtkSrc in tvtkGridFunctionSrcs]
            except TypeError:
                self.tvtkGridFunctionSrcs = [VTKDataSource(data=tvtkGridFunctionSrcs)]
        else:
            self.tvtkGridFunctionSrcs = []
        if tvtkStructuredGridDataSrcs is not None:
            try:
                self.tvtkStructuredGridDataSrcs = [VTKDataSource(data=tvtkSrc)
                                             for tvtkSrc in tvtkStructuredGridDataSrcs]
            except TypeError:
                self.tvtkStructuredGridDataSrcs = [VTKDataSource(data=tvtkStructuredGridDataSrcs)]
        else:
            self.tvtkStructuredGridDataSrcs = []
        self.dataRange = dataRange

    @on_trait_change('scene.activated')
    def create_plot(self):
        from mayavi.modules.api import Surface, Vectors
        self.engine = self.scene.engine
        self.surface = []
        self.surface1 = []
        # for src in self.tvtkGridFunctionSrcs:
        #     self.surface.append(Surface())
        #     self.surface1.append(Surface())
        #     self.surface1[-1].actor.property.representation = 'wireframe'
        #     self.surface1[-1].actor.actor.visibility = self.enable_grid
        #     self.surface[-1].actor.actor.visibility=self.enable_surface
        #     self.engine.add_source(src)
        #     self.engine.add_module(self.surface[-1], obj=src)
        #     self.engine.add_module(self.surface1[-1], obj=src)
        #     self.surface[-1].actor.actor.visibility=self.enable_surface
        for src in self.tvtkGridFunctionSrcs + self.tvtkStructuredGridDataSrcs:
            self.surface.append(Surface())
            self.surface[-1].actor.actor.visibility=self.enable_surface
            self.engine.add_source(src)
            self.engine.add_module(self.surface[-1], obj=src)
        for src in self.tvtkGridSrcs:
            self.surface1.append(Surface())
            self.surface1[-1].actor.property.representation = 'wireframe'
            self.surface1[-1].actor.actor.visibility = self.enable_grid
            self.surface1[-1].actor.mapper.scalar_visibility = False
            self.engine.add_source(src)
            self.engine.add_module(self.surface1[-1], obj=src)
        self.module_manager = self.engine.scenes[0].children[0].children[0]
        if self.dataRange is not None:
            self.module_manager.scalar_lut_manager.data_range = self.dataRange
            self.module_manager.scalar_lut_manager.use_default_range = False
        self.module_manager.scalar_lut_manager.show_legend = self.enable_legend
        if self.point_cell == "Point Data":
            self.module_manager.lut_data_mode = 'point data'
        else:
            self.module_manager.lut_data_mode = 'cell data'

    @on_trait_change('real_imag')
    def update_real_imag(self):
        if self.real_imag=="Real Part":
            for src in self.tvtkGridFunctionSrcs + self.tvtkStructuredGridDataSrcs:
                src.point_scalars_name = 'real'
                try:
                    src.cell_scalars_name = 'real'
                except TraitError:
                    pass # cell data not present
        elif self.real_imag=="Imaginary Part":
            for src in self.tvtkGridFunctionSrcs + self.tvtkStructuredGridDataSrcs:
                src.point_scalars_name = 'imag'
                try:
                    src.cell_scalars_name = 'imag'
                except TraitError:
                    pass # cell data not present
        else:
            for src in self.tvtkGridFunctionSrcs + self.tvtkStructuredGridDataSrcs:
                src.point_scalars_name = 'abs^2'
                try:
                    src.cell_scalars_name = 'abs^2'
                except TraitError:
                    pass # cell data not present

    @on_trait_change('enable_grid')
    def update_grid(self):
        for s in self.surface1:
            s.actor.actor.visibility = self.enable_grid

    @on_trait_change('point_cell')
    def update_point_cell(self):
        if self.point_cell == "Point Data":
            self.module_manager.lut_data_mode = 'point data'
        else:
            self.module_manager.lut_data_mode = 'cell data'

    @on_trait_change('enable_surface')
    def update_surface(self):
        for s in self.surface:
            s.actor.actor.visibility = self.enable_surface

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
                    orientation="horizontal")),
                orientation="vertical"),
                resizable=True,title="Grid Function Viewer")

def tvtkGrid(grid):
    """Return a TVTK object visualizing the specified grid."""

    if grid.topology()=="triangular":
        (points,elems,auxData) = grid.leafView().getRawElementData()
        elem_list = elems[:-1,:].T
        mesh = tvtk.PolyData()
        mesh.points = points.T
        mesh.polys = elem_list
    else:
        raise TypeError("Visualization of this grid topology not implemented!")
    return mesh

def tvtkGridFunction(g):
    """Return a TVTK object visualizing the specified grid function."""
    tvtkObj = tvtkGrid(g.grid())
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
        tvtkObj.point_data.set_active_vectors('real')
    elif g.componentCount()==1:
        tvtkObj.cell_data.set_active_scalars('real')
        tvtkObj.point_data.set_active_scalars('real')
    else:
        raise Exception("plotGridFunction: Only GridFunctions with "
                        "componentCount 1 or 3 are supported.")

    return tvtkObj

def tvtkStructuredGridData(points, data, dims):
    """Return a TVTK object visualizing data on a structured grid.

    *Parameters:*
        - points (2D NumPy array):
            Array of point coordinates of dimension (3, np), where np is the
            number of points. The points should be lying on a regular grid.
        - data (1D or 2D NumPy array)
            Array of data values corresponding to each point. It can be either a
            1D array of length np (number of points) or a 2D array of dimensions
            (1, np) or (3, np).
        - dims
            A tuple of two numbers representing the first and second dimension
            of the grid of points. Their product must be equal to np.

    Example usage::

        import numpy as np
        from bempp import visualization2 as vis

        evaluationOptions = createEvaluationOptions()

        x, y, z = np.mgrid[0:1:5j, 0:2:10j, 0:0:1j]
        nx, ny = x.shape[:2]
        points = np.vstack((x.ravel(), y.ravel(), z.ravel()))
        field = somePotentialOperator.evaluateAtPoints(someGridFunction, points,
                                                       evaluationOptions)
        tvtkObj1 = vis.tvtkStructuredGridData(points, field, (nx, ny))

        x, y, z = np.mgrid[0:0:1j, 0:2:10j, 0:5:7j]
        ny, nz = x.shape[1:]
        points = np.vstack((x.ravel(), y.ravel(), z.ravel()))
        field = somePotentialOperator.evaluateAtPoints(someGridFunction, points,
                                                       evaluationOptions)
        tvtkObj2 = vis.tvtkStructuredGridData(points, field, (ny, nz))
        plotScalarData(tvtkStructuredData=[tvtkObj1, tvtkObj2])
    """

    points = np.asanyarray(points)
    data = np.asanyarray(data)
    if points.ndim != 2:
        raise ValueError("'points' must be a 2D array")
    dimWorld, pointCount = points.shape
    if dimWorld != 3:
        raise ValueError("'points' must have 3 coordinates")
    if data.ndim == 1:
        data = data[np.newaxis,:]
    if data.ndim != 2:
        raise ValueError("'data' must be a 1D or 2D array")
    componentCount = data.shape[0]
    if componentCount not in (1, 3):
        raise ValueError("'data' must have either 1 or three rows")
    if data.shape[1] != pointCount:
        raise ValueError("The numbers of points in 'points' and 'data'"
                         " do not match")
    if len(dims) != 2:
        raise ValueError("'dims' must be a two-element tuple")
    if dims[0] * dims[1] != pointCount:
        raise ValueError("The product of dims[0] and dims[1] does not match "
                         "the point count")
    tvtkObj = tvtk.StructuredGrid(dimensions=dims[::-1] + (1,))
    tvtkObj.points = points.T
    tvtkObj.point_data.add_array(np.real(data.T))
    tvtkObj.point_data.add_array(np.imag(data.T))
    tvtkObj.point_data.add_array(np.sum(abs(data)**2,axis=0))
    tvtkObj.point_data.get_abstract_array(0).name = 'real'
    tvtkObj.point_data.get_abstract_array(1).name = 'imag'
    tvtkObj.point_data.get_abstract_array(2).name = 'abs^2'

    if componentCount == 3:
        tvtkObj.point_data.set_active_scalars('abs^2')
        tvtkObj.point_data.set_active_vectors('real')
    else:
        tvtkObj.point_data.set_active_scalars('real')
    return tvtkObj

def plotGrid(grid, representation='wireframe'):
    """Visualize a grid.

    *Parameters:*
        - gf (Grid):
            The grid to be plotted.
        - representation (str)
            Plot type (e.g. 'wireframe', 'surface').

    Returns a Traits object that contains the visualization."""

    gridTvtkData = tvtkGrid(grid)
    return mlab.pipeline.surface(gridTvtkData, representation=representation)

def plotGridFunction(gf, dataRange=None):
    """Visualize a grid function.

    *Parameters:*
        - gf (GridFunction):
             The grid function to be plotted.

    Returns a Traits object that contains the visualization."""

    tvtkObj = tvtkGridFunction(gf)
    if gf.componentCount()==3:
        v = _VectorVisualization(tvtkObj, tvtkObj, None, dataRange)
    else:
        v = _ScalarVisualization(tvtkObj, tvtkObj, None, dataRange)
    v.configure_traits()
    return v

def plotStructuredGridData(points, data, dims, dataRange=None):
    """Visualize data sampled on a regular 2D grid of points.

    *Parameters:*
        - points (2D NumPy array):
            Array of point coordinates of dimension (3, np), where np is the
            number of points. The points should be lying on a regular grid.
        - data (1D or 2D NumPy array)
            Array of data values corresponding to each point. It can be either a
            1D array of length np (number of points) or a 2D array of
            dimensions (1, np) or (3, np).
        - dims
           A tuple of two numbers representing the first and second dimension
           of the grid of points. Their product must be equal to np.
        - dataRange
           Can be None or a tuple of two floats. In the first case, the data
           range of the plot will be determined automatically, otherwise
           it will be set to the specified tuple.
            
    Returns a Traits object that contains the visualization.

    Example usage::

        import numpy as np
        from bempp import visualization2 as vis

        x, y, z = np.mgrid[0:1:5j, 0:2:10j, 0:0:1j]
        nx, ny = x.shape[:2]
        points = np.vstack((x.ravel(), y.ravel(), z.ravel()))

        evaluationOptions = createEvaluationOptions()
        field = somePotentialOperator.evaluateAtPoints(someGridFunction, points,
                                                       evaluationOptions)
        vis.plotStructuredGridData(points, field, (nx, ny))
    """

    tvtkObj = tvtkStructuredGridData(points, data, dims)
    if data.ndim == 2 and data.shape[0] == 3:
        v = _VectorVisualization(None, None, tvtkObj, dataRange)
    else:
        v = _ScalarVisualization(None, None, tvtkObj, dataRange)
    v.configure_traits()
    return v

def plotScalarData(tvtkGrids=None, tvtkGridFunctions=None, tvtkStructuredGridData=None,
                   dataRange=None):
    """Visualize multiple scalar data sets in the same plot.

    This function plots an arbitrary number of grids, scalar-valued grid functions and
     data sets sampled on regular 2D grids of points in a single window.

    *Parameters:*
        - tvtkGrids:
            None, an object returned by tvtkGrid(), or a list of such
            objects.
        - tvtkGridFunctions:
            None, an object returned by tvtkGridFunction(), or a list of
            such objects.
        - tvtkStructuredGridData
            None, an object returned by tvtkStructuredGridData(), or a list
            of such objects.
        - dataRange
            Can be None or a tuple of two floats. In the first case, the data
            range of the plot will be determined automatically, otherwise
            it will be set to the specified tuple.

    Returns a Traits object that contains the visualization.
    """
    return plotData(tvtkGrids, tvtkGridFunctions, tvtkStructuredGridData,
                    dataRange, scalar=True)

def plotVectorData(tvtkGrids=None, tvtkGridFunctions=None, tvtkStructuredGridData=None,
                   dataRange=None):
    """Visualize multiple scalar data sets in the same plot.

    This function plots an arbitrary number of grids, vector-valued grid functions and
    data sets sampled on regular 2D grids of points in a single window.

    *Parameters:*
        - tvtkGrids:
            None, an object returned by tvtkGrid(), or a list of such
            objects.
        - tvtkGridFunctions:
            None, an object returned by tvtkGridFunction(), or a list of
            such objects.
        - tvtkStructuredGridData
            None, an object returned by tvtkStructuredGridData(), or a list
            of such objects.
        - dataRange
            Can be None or a tuple of two floats. In the first case, the data
            range of the plot will be determined automatically, otherwise
            it will be set to the specified tuple.

    Returns a Traits object that contains the visualization.
    """
    return plotData(tvtkGrids, tvtkGridFunctions, tvtkStructuredGridData,
                    dataRange, scalar=False)

def plotData(tvtkGrids=None, tvtkGridFunctions=None, tvtkStructuredGridData=None,
             dataRange=None, scalar=True):
    """Visualize multiple scalar or vector data sets in the same plot.

    This function plots an arbitrary number of grids, grid functions and data
    sets sampled on regular 2D grids of points in a single window. The grid
    functions and structured grid data can be scalar- or vector-valued depending
    on the value of the parameter 'scalar'.

    *Parameters:*
        - tvtkGrids:
            None, an object returned by tvtkGrid(), or a list of such
            objects.
        - tvtkGridFunctions:
            None, an object returned by tvtkGridFunction(), or a list of
            such objects.
        - tvtkStructuredGridData
            None, an object returned by tvtkStructuredGridData(), or a list
            of such objects.
        - dataRange
            Can be None or a tuple of two floats. In the first case, the data
            range of the plot will be determined automatically, otherwise
            it will be set to the specified tuple.
        - scalar (bool)
            If True, the arguments tvtkGridFunctions and tvtkStructuredGridData
            should represent scalar-valued data sets, otherwise vector-valued
            data sets.

    Returns a Traits object that contains the visualization.
    """
    if tvtkGrids is None:
        tvtkGrids = tvtkGridFunctions
    if scalar:
        v = _ScalarVisualization(tvtkGrids, tvtkGridFunctions, tvtkStructuredGridData,
                                 dataRange)
    else:
        v = _VectorVisualization(tvtkGrids, tvtkGridFunctions, tvtkStructuredGridData,
                                 dataRange)
    v.configure_traits()
    return v

def plotThreePlanes(points, vals):
    """Plot a three planes view of data.

       *Parameters:*
           - points:
               An mgrid object defining points in a box.
           - vals:
               The corresponding values.

       Potential values in the correct format are returned by the function
       tools.evaluatePotentialInBox.
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

