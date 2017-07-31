"""Define routines specific to IPython and IPython notebooks."""

import numpy as _np

class IPyNotebookSurfaceViewer(object):
    """A surface viewer for IPython Notebooks"""

    def __init__(self, vertices, elements, data=None):
        """
        Initialize the viewer

        Attributes
        ----------
        vertices : np.ndarray
            A (3xN) array of coordinates for N vertices
        elements : np.ndarray
            A (3xM) array defining M elements
        data : np.ndarray
            A (M x 3 x dim) array of data. data[i, j, k] is the kth component
            of the data associated with the jth node for the ith face.

        """
        self._vertices = vertices.T
        self._elements = elements.T

        self._data = data
        self._mesh = None
        self._wireframe = None
        self._scene = None
        self._renderer = None
        self._cmap = None
        self._vis_data = None
        self._vmin = None
        self._vmax = None
                
    @property
    def data(self):
        """Return the data."""
        return self._data

    @property
    def is_scalar(self):
        """Return true if data is scalar."""
        return self.data.shape[2] == 1
                        
    @property
    def is_real(self):
        """Return true if data is real."""
        return _np.all(_np.isreal(self.data))

    @property
    def vertices(self):
        """Return the vertices."""
        return self._vertices

    @property
    def elements(self):
        """Return the elements."""
        return self._elements

    def visualize(self):
        """Start the visualization and initialize widgets"""
        from pythreejs import PlainGeometry, Mesh, LambertMaterial, \
                PhongMaterial, DirectionalLight, \
                PerspectiveCamera, Scene, AmbientLight, \
                Renderer, OrbitControls
        from matplotlib import colors
        from matplotlib import cm
        from IPython.display import display
        import numpy as np
        import ipywidgets as widgets


        bbox = self._compute_bounding_box()

        center = .5 * _np.sum(bbox, axis=0)
        diam =  bbox[1, :] - bbox[0, :]

        position = (center + 2.5 * diam).tolist()

        if self.data is not None:
            vis_data = self.data
            if not self.is_real:
                vis_data = _np.real(self.data)
            if not self.is_scalar:
                vis_data = _np.sum(_np.abs(self.data)**2, axis=2)
            self._vis_data = vis_data
            self._cmap = cm.jet
            self._vmin = vis_data.min()
            self._vmax = vis_data.max()
            cnorm = colors.Normalize(vmin=self._vmin, vmax=self._vmax)
            facecolors = self._convert_data_to_colors(self._cmap, cnorm, vis_data)
        else:
            facecolors = self.elements.shape[0] * [3 * [_np.array([1., 1., 1.])]]
        self._geom = PlainGeometry(vertices=self.vertices, faces=self.elements, faceColors = facecolors)

        self._mesh = Mesh(geometry=self._geom, material = LambertMaterial(vertexColors = 'VertexColors'))
        self._wireframe = Mesh(geometry=self._geom, material = PhongMaterial(wireframe=True, color='black'))
        light = DirectionalLight(color='white', position=position, intensity=0.5)
        camera = PerspectiveCamera(position=position, fov=20)
        self._scene = Scene(children=[self._mesh, self._wireframe, AmbientLight(color='white')])

        self._renderer = Renderer(camera=camera, background='white', background_opacity=1,
                                scene = self._scene, controls=[OrbitControls(controlling=camera)])

        if self.data is not None:
 
            # Enable/Disable wireframe
            wireframe_toggle = widgets.Checkbox(
            value=self._wireframe.visible,
            description='Enable wireframe',
            disabled=False
            )
            wireframe_toggle.observe(self.on_toggle_wireframe, names='value')

            # Change vmin/vmax
            vmin_box = widgets.FloatText(
                value=self._vmin,
                description='vmin:',
                disabled=False
            ) 
            vmin_box.observe(self.on_change_vmin, names='value')

            vmax_box = widgets.FloatText(
                value=self._vmax,
                description='vmax:',
                disabled=False
            ) 
            vmax_box.observe(self.on_change_vmax, names='value')

            vmin_info = widgets.Label(
                    value='Lower bound: {0}'.format(self._vmin))
            vmax_info = widgets.Label(
                    value='Upper bound: {0}'.format(self._vmax))

            range_info_box = widgets.VBox([vmin_info, vmax_info])
            range_change = widgets.HBox([vmin_box, vmax_box])
            vbox = widgets.VBox([range_info_box, range_change, wireframe_toggle])
            display(self._renderer, vbox)       
        else:
            display(self._renderer) 


    def update_geometry_color(self):
        """Update color of geometry"""
        from matplotlib import colors
        cnorm = colors.Normalize(vmin=self._vmin, vmax=self._vmax)
        self._geom.faceColors = self._convert_data_to_colors(
                self._cmap, cnorm, self._vis_data)

    def on_change_vmin(self, change):
        """Toggle vmin"""
        if change.new > self._vmax:
            change.owner.value = change.old
            return
        self._vmin = change.new
        self.update_geometry_color()

    def on_change_vmax(self, change):
        """Toggle vmin"""
        if change.new < self._vmin:
            change.owner.value = change.old
            return
        self._vmax = change.new
        self.update_geometry_color()


    def on_toggle_wireframe(self, change):
        """Toggle wireframe change."""
        self._wireframe.visible = change.new

    def _convert_data_to_colors(self, cmap, cnorm, data):
        return [[cmap(cnorm(value))[:3] for value in elem.flat] for elem in self.data]


    def _compute_bounding_box(self):

        xmin = self.vertices[:, 0].min()
        xmax = self.vertices[:, 0].max()
        ymin = self.vertices[:, 1].min()
        ymax = self.vertices[:, 1].max()
        zmin = self.vertices[:, 2].min()
        zmax = self.vertices[:, 2].max()

        return _np.array([[xmin, ymin, zmin], [xmax, ymax, zmax]])









