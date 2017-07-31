"""Define interfaces to external viewers."""
import numpy as _np

def visualize(obj):
    """Main visualization method"""
    import bempp.api
    
    if bempp.api.PLOT_BACKEND == 'gmsh':
        visualize_with_gmsh(obj)
    if bempp.api.PLOT_BACKEND == 'ipython_notebook':
        visualize_with_ipython_notebook(obj)

def visualize_with_gmsh(obj):
    """
    View a grid or grid function with Gmsh

    Parameters
    ----------
    obj : bempp.api.Grid or bempp.api.GridFunction
        Grid or grid function to visualize.

    Notes
    -----
    This function writes the data into a temp file and
    visualizes it.

    """
    import tempfile
    import subprocess
    from bempp.api import export, GMSH_PATH, TMP_PATH, GridFunction
    from bempp.api.grid.grid import Grid
    from numpy import real

    if GMSH_PATH is None:
        print("Gmsh not available for visualization.")
        return None

    outfile = tempfile.NamedTemporaryFile(
        suffix=".msh", dir=TMP_PATH, delete=False)
    if isinstance(obj, Grid):
        export(grid=obj, file_name=outfile.name)
    elif isinstance(obj, GridFunction):
        export(grid_function=obj, file_name=outfile.name, transformation=real)
    outfile.close()

    subprocess.Popen([GMSH_PATH, outfile.name])

def visualize_with_ipython_notebook(obj):
    """View a grid or grid function in an IPython Notebook"""
    from bempp.api.utils.ipython import IPyNotebookSurfaceViewer
    from bempp.api import GridFunction
    from bempp.api.grid.grid import Grid

    if isinstance(obj, Grid):
        viewer = IPyNotebookSurfaceViewer(obj.leaf_view.vertices,
            obj.leaf_view.elements)
        viewer.visualize()
    elif isinstance(obj, GridFunction):
        grid = obj.space.grid
        index_set = grid.leaf_view.index_set()
        local_coordinates = _np.array([[0, 1, 0], [0, 0, 1]])
        values = _np.zeros((grid.leaf_view.entity_count(0), 3, 
                            obj.component_count), dtype=obj.dtype)
        for index, element in enumerate(grid.leaf_view.entity_iterator(0)):
            local_values = obj.evaluate(element, local_coordinates)
            values[index,:,:] = local_values.T
        viewer = IPyNotebookSurfaceViewer(grid.leaf_view.vertices,
            grid.leaf_view.elements, values)
        viewer.visualize()

def set_gmsh_viewer():
    """Change plotting default to Gmsh."""
    import bempp.api
    bempp.api.PLOT_BACKEND = 'gmsh'

def set_ipython_notebook_viewer():
    """Change plotting default to IPython"""
    import bempp.api
    bempp.api.PLOT_BACKEND = 'ipython_notebook'
       
