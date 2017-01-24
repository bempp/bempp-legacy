from structure cimport init_structure, dealloc_structure
from cpython.mem cimport PyMem_Malloc, PyMem_Free

cdef class Structure:
    """ A useless class """
    def __cinit__(self):
        self.cdata = <CStructure*> PyMem_Malloc(sizeof(CStructure))
        if not self.cdata:
            raise MemoryError("Could not allocate C object")
        init_structure(self.cdata)
    def __dealloc__(self):
        dealloc_structure(self.cdata)
        PyMem_Free(self.cdata)

    property meaning_of_life:
        def __get__(self):
            return self.cdata.meaning_of_life
        def __set__(self, value):
            self.cdata.meaning_of_life = int(value)
    property message:
        def __get__(self):
            return self.cdata.message
