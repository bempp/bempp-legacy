

cdef class CudaGrid:

    def __cinit__(self):
        self.impl_.reset()

    def __init__(self):
        pass

    def __dealloc__(self):
        self.impl_.reset()


