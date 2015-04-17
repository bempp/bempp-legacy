from cython.operator cimport dereference as deref
from cython.operator cimport address
from libcpp.string cimport string 
from libcpp cimport bool as cbool
from .byte_conversion import convert_to_bytes
from libcpp.vector cimport vector

cdef class ParameterList:

    def __cinit__(self):
        self.impl_ = new c_ParameterList()

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.impl_

    




        

        



