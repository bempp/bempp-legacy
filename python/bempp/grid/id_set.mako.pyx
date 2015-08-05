from cython.operator cimport dereference as deref
from bempp.grid.entity cimport Entity0, Entity1, Entity2

cdef class IdSet:

    def __cinit__(self):
        pass

    def __init__(self):
        pass

    def __dealloc__(self):
        pass

    cdef size_t entity_id_0(self, Entity0 entity):
        return deref(self.impl_).entityId(deref(entity.impl_))

    cdef size_t entity_id_1(self, Entity1 entity):
        return deref(self.impl_).entityId(deref(entity.impl_))

    cdef size_t entity_id_2(self, Entity2 entity):
        return deref(self.impl_).entityId(deref(entity.impl_))

    def entity_id(self, entity):

        if isinstance(entity, Entity0):
            return self.entity_id_0(entity)
        elif isinstance(entity, Entity1):
            return self.entity_id_1(entity)
        elif isinstance(entity, Entity2):
            return self.entity_id_2(entity)
        else:
            raise ValueError("Unknown type {0}".format(str(type(entity))))

    cpdef size_t sub_entity_id(self, Entity0 element, size_t i, int codim_sub):
        return deref(self.impl_).subEntityId(deref(element.impl_), i, codim_sub)






