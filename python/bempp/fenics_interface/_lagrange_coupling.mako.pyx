from bempp.space.space cimport Space
from cython.operator cimport dereference as deref
cimport numpy as np


def p1_vertex_map(Space space):

    return _py_p1_vertex_map(deref(space.impl_))


