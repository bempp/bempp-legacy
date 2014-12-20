from libcpp import bool as cbool
from cython.operator import dereference as deref

cdef extern from "<memory>" namespace "std":
    cdef cppclass unique_ptr[T]:
        unique_ptr()
        unique_ptr(T*)
        T& operator*()
        T* get()
        T* release()
        void reset(T*)
        void swap(unique_ptr[T]& other)
