from bempp.core.fmm.octree cimport Octree
from bempp.core.utils.shared_ptr cimport shared_ptr
from cython.operator cimport dereference as deref



cdef class InteractionList:

    def __cinit__(self, Octree octree, unsigned long node_index, unsigned int level):
        pass

    def __init__(self, Octree octree, unsigned long node_index, unsigned int level):
        self.impl_.assign(shared_ptr[c_InteractionList](new c_InteractionList(
            deref(octree.impl_), node_index, level)))

    def __dealloc__(self):
        self.impl_.reset()

    def __iter__(self):
        return self

    def next(self):
        """Move to next element."""
        if not self._finished():
            val =  deref(deref(self.impl_))
            deref(self.impl_).next()
            return val
        else:
            raise StopIteration()

    def __next__(self):
        return self.next()

    def _finished(self):
        """Check if finished."""
        return deref(self.impl_).finished()


