"""Define interfaces to external viewers."""

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
