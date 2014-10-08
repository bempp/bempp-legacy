""" Grids """
__all__ = ['Grid', 'sample']
from .grid import Grid


def sample(topology='triangular', filename="sphere-h-0.4.msh"):
    """ Convenience function to access the sample grids """
    from os.path import join, exists
    from ..config import paths
    filename = join(paths.meshes, filename)
    if not exists(filename):
        raise IOError("Mesh %s does not exist" % filename)
    return Grid(topology="triangular", filename=filename)
