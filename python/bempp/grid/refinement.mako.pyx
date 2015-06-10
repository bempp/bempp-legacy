from .grid cimport Grid, c_Grid
from .entity cimport Entity0
from bempp.utils cimport shared_ptr, const_pointer_cast, reverse_const_pointer_cast
from cython.operator cimport dereference as deref

cdef class RefinementFactory:
    cdef shared_ptr[c_Grid] grid_ptr

    def __cinit__(self,Grid grid):
        pass

    def __dealloc__(self):
        self.grid_ptr.reset()

    def __init__(self,Grid grid):
        self.grid_ptr = const_pointer_cast[c_Grid]((<Grid>grid.clone()).impl_)
        

    def mark(self, Entity0 element, modus='refine'):
        cdef int ref_count

        if modus == 'refine':
            ref_count = 1
        elif modus == 'coarsen':
            ref_count = -1
        return deref(self.grid_ptr).mark(ref_count, deref(element.impl_))

    def get_mark(self, Entity0 element):

        return deref(self.grid_ptr).getMark(deref(element.impl_))

    def refine(self):
        
        deref(self.grid_ptr).preAdapt()
        deref(self.grid_ptr).adapt()
        deref(self.grid_ptr).postAdapt()
        cdef Grid grid = Grid.__new__(Grid)
        grid.impl_ = reverse_const_pointer_cast[c_Grid](self.grid_ptr)
        return grid

    def global_refine(self, int ref_count):

        deref(self.grid_ptr).globalRefine(ref_count)
        cdef Grid grid = Grid.__new__(Grid)
        grid.impl_ = reverse_const_pointer_cast[c_Grid](self.grid_ptr)
        return grid
        

    
        
