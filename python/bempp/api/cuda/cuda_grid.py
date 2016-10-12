

class CudaGrid(object):
    """A class for grids living on a CUDA capable device."""
    
    def __init__(self, impl):
        self._impl = impl
        
    def setup_geometry(self):
        """Setup the geometry on the device."""
        self._impl.setup_geometry()

