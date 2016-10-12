from cython.operator cimport dereference as deref

cdef class CudaGrid:
    """A class for grids living on a CUDA capable device."""
    
    def __cinit__(self):
        self.impl_.reset()

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()
        
    def setup_geometry(self):
        """ Setup geometry on the device."""
        deref(self.impl_).setupGeometry()


