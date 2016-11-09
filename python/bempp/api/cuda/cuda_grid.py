

class CudaGrid(object):
    """A class for grids living on a CUDA capable device."""
    
    def __init__(self, impl):
        self._impl = impl

