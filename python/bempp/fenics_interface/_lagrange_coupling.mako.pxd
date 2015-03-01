
from bempp.utils cimport shared_ptr
from bempp.space.space cimport c_Space
from bempp.space.space cimport SpaceVariants

cdef extern from "bempp/fenics_interface/py_lagrange_interface.hpp" namespace "Bempp":
    object _py_p1_vertex_map(const SpaceVariants)

