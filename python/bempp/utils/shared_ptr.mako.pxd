from libcpp import bool as cbool
from cython.operator import dereference as deref

cdef extern from "<boost/shared_ptr.hpp>" namespace "boost":
    cdef cppclass shared_ptr[T]:
        shared_ptr()
        shared_ptr(T*)
        shared_ptr(const shared_ptr[T]&)
        void assign "operator="(shared_ptr[T]) 
        T& operator*()
        T* get()
        void reset(T*) except *

    shared_ptr[T] static_pointer_cast[T,U](const shared_ptr[U]&)
