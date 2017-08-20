"""Define interfaces to external viewers."""
import numpy as _np

def visualize(obj, mode='vertices', transformation=None):
    """
    Main visualization method

    Attributes
    ----------
    obj : Grid or GridFunction object
        A grid or grid function to plot
    mode : string
        One of 'element' or 'node'. If 'element' is chosen
        the color is determined by the mid-point of the faces
        of the grid. For 'node' the vertex values are
        chosen (default: 'element'). Only used for
        grid functions.
    transformation : string or object
        One of 'real', 'imag', 'abs', 'log_abs' or
        'abs_squared' or a callable object. 
        Describes the data transformation
        before plotting. For functions with vector values
        only 'abs', 'log_abs' or 'abs_squared' are allowed.
        If a callable object is given this is applied instead.
        It is important that the callable returns numpy arrays
        with the same number of dimensions as before.

    """
    import bempp.api
    import numpy as np

    transform = None

    if transformation == 'real':
        transform = np.real
    elif transformation == 'imag':
        transform = np.imag
    elif transformation == 'abs':
        transform = lambda x: np.sqrt(np.sum(np.abs(x)**2, axis=0, keepdims=True))
    elif transformation == 'log_abs':
        transform = lambda x: np.log(np.sqrt(
            np.sum(np.abs(x)**2, axis=0, keepdims=True)))
    elif transformation == 'abs_squared':
        transform = lambda x: np.sum(np.abs(x)**2, axis=0, keepdims=True)
    else:
        transform = transformation

    if bempp.api.PLOT_BACKEND == 'gmsh':
        visualize_with_gmsh(obj, mode, transform)
    if bempp.api.PLOT_BACKEND == 'paraview':
        visualize_with_paraview(obj, mode, transform)
    if bempp.api.PLOT_BACKEND == 'ipython_notebook':
        visualize_with_ipython_notebook(obj, mode, transform)

def visualize_with_gmsh(obj, mode='element', transformation=None):
    """
    View a grid or grid function with Gmsh

    Parameters
    ----------
    obj : bempp.api.Grid or bempp.api.GridFunction
        Grid or grid function to visualize.
    mode : string
        One of 'element' or 'node'
        (default 'vertices')
    transformation : callable
        A function object that is applied to the data before
        writing it out

    Notes
    -----
    This function writes the data into a temp file and
    visualizes it.

    """
    import tempfile
    import subprocess
    from bempp.api import export, GMSH_PATH, TMP_PATH, GridFunction
    from bempp.api.grid.grid import Grid
    import numpy as np

    if transformation is None:
        transformation = np.real

    if GMSH_PATH is None:
        print("Gmsh not available for visualization.")
        return None


    outfile = tempfile.NamedTemporaryFile(
        suffix=".msh", dir=TMP_PATH, delete=False)
    if isinstance(obj, Grid):
        export(grid=obj, file_name=outfile.name)
    elif isinstance(obj, GridFunction):
        export(grid_function=obj, file_name=outfile.name, 
                transformation=transformation, data_type=mode)
    outfile.close()

    subprocess.Popen([GMSH_PATH, outfile.name])

def visualize_with_paraview(obj, mode='element', transformation=None):
    """
    View a grid or grid function with Paraview

    Parameters
    ----------
    obj : bempp.api.Grid or bempp.api.GridFunction
        Grid or grid function to visualize.
    mode : string
        One of 'element' or 'node'
        (default 'vertices')
    transformation : callable
        A function object that is applied to the data before
        writing it out

    Notes
    -----
    This function writes the data into a temp file and
    visualizes it.

    """
    import tempfile
    import subprocess
    from bempp.api import export, GMSH_PATH, TMP_PATH, GridFunction
    from bempp.api.grid.grid import Grid
    import numpy as np
    from bempp.api.utils import which

    pview = which("paraview")
    if pview is None:
        raise EnvironmentError(
                "Could not find Paraview." +
                "Interactive plotting with Paraview not available.")

    if transformation is None:
        transformation = np.real

    outfile = tempfile.NamedTemporaryFile(
        suffix=".vtk", dir=TMP_PATH, delete=False)
    print(outfile)
    if isinstance(obj, Grid):
        export(grid=obj, file_name=outfile.name)
    elif isinstance(obj, GridFunction):
        export(grid_function=obj, file_name=outfile.name, 
                transformation=transformation, data_type=mode)
    outfile.close()

    subprocess.Popen([pview, outfile.name])

def visualize_with_ipython_notebook(obj, mode='element', transformation=None):
    """View a grid or grid function in an IPython Notebook"""
    from bempp.api.utils.ipython import IPyNotebookSurfaceViewer
    from bempp.api import GridFunction
    from bempp.api.grid.grid import Grid
    import numpy as np

    if transformation is None:
        transformation = np.real

    if isinstance(obj, Grid):
        viewer = IPyNotebookSurfaceViewer(obj.leaf_view.vertices,
            obj.leaf_view.elements)
        viewer.visualize()
    elif isinstance(obj, GridFunction):
        grid = obj.space.grid
        index_set = grid.leaf_view.index_set()
        if mode == 'element':
            local_coordinates = _np.array([[0, 1, 0], [0, 0, 1]])
        else:
            local_coordinates = _np.array([[1./3, 1./3, 1./3],
                [1./3, 1./3, 1./3]])
        values = _np.zeros((grid.leaf_view.entity_count(0), 3, 
                            obj.component_count), dtype=obj.dtype)
        for element in grid.leaf_view.entity_iterator(0):
            index = index_set.entity_index(element)
            local_values = transformation(obj.evaluate(element, local_coordinates))
            values[index,:,:] = local_values.T
        viewer = IPyNotebookSurfaceViewer(grid.leaf_view.vertices,
            grid.leaf_view.elements, values)
        viewer.visualize()

def set_gmsh_viewer():
    """Change plotting default to Gmsh."""
    import bempp.api
    bempp.api.PLOT_BACKEND = 'gmsh'

def set_paraview_viewer():
    """Change plotting default to Paraview."""
    import bempp.api
    bempp.api.PLOT_BACKEND = 'paraview'

def set_ipython_notebook_viewer():
    """Change plotting default to IPython"""
    import bempp.api
    bempp.api.PLOT_BACKEND = 'ipython_notebook'
       
