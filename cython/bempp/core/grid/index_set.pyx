from cython.operator cimport dereference as deref
from bempp.core.grid.entity cimport Entity0, Entity1, Entity2

cdef class IndexSet:

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        pass

    cdef size_t entity_index_0(self, Entity0 entity):
        return deref(self.impl_).entityIndex(deref(entity.impl_))

    cdef size_t entity_index_1(self, Entity1 entity):
        return deref(self.impl_).entityIndex(deref(entity.impl_))

    cdef size_t entity_index_2(self, Entity2 entity):
        return deref(self.impl_).entityIndex(deref(entity.impl_))

    def entity_index(self, entity):

        if isinstance(entity, Entity0):
            return self.entity_index_0(entity)
        elif isinstance(entity, Entity1):
            return self.entity_index_1(entity)
        elif isinstance(entity, Entity2):
            return self.entity_index_2(entity)
        else:
            raise ValueError("Unknown type {0}".format(str(type(entity))))

    cpdef size_t sub_entity_index(self, Entity0 element, size_t i, int codim_sub):
        return deref(self.impl_).subEntityIndex(deref(element.impl_), i, codim_sub)






