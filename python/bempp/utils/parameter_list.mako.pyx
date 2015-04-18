from cython.operator cimport dereference as deref
from cython.operator cimport address
from libcpp.string cimport string 
from libcpp cimport bool as cbool
from .byte_conversion import convert_to_bytes
from libcpp.vector cimport vector

cdef class _AssemblyParameterList:

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        # Pointer deallocated by parent class
        pass

    property boundary_assembly_type:

        def __get__(self):

            cdef char* s = b"options.assembly.boundaryOperatorAssemblyType"
            return (deref(self.impl_).get_string(s)).decode("UTF-8")

        def __set__(self,object value):

            cdef char* s = b"options.assembly.boundaryOperatorAssemblyType"           
            cdef string stringVal = convert_to_bytes(value)

            deref(self.impl_).put_string(s,stringVal)



cdef class ParameterList:

    def __cinit__(self):
        self.impl_ = new c_ParameterList()
        self._assembly = _AssemblyParameterList()
        (<_AssemblyParameterList>self._assembly).impl_ = self.impl_

    def __init__(self):
        pass

    def __dealloc__(self):
        del self.impl_

    property assembly:

        def __get__(self):
            return self._assembly





        

        



