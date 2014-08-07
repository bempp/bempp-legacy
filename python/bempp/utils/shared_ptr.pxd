from cython.operator import dereference as deref

cdef extern from "<boost/shared_ptr.hpp>" namespace "boost":
    cdef cppclass shared_ptr[T]:
        T* get()
        T& operator*()
