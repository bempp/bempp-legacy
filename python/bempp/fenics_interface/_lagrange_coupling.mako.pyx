from bempp.space.space cimport Space
cimport numpy as np

np.import_array()

def p1_vertex_map(Space space):

    return _py_p1_vertex_map(space.impl_)


