from cython.operator cimport dereference as deref

cdef class CudaGridSingle:
    """A class for grids living on a CUDA capable device."""
   
    def __cinit__(self):
        self.impl_.reset()

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()


cdef class CudaGridDouble:
    """A class for grids living on a CUDA capable device."""
   
    def __cinit__(self):
        self.impl_.reset()

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()
