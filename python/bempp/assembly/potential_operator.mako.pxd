from bempp.assembly.discrete_boundary_operator cimport DiscreteBoundaryOperatorBase
from bempp.space.space cimport Space
import numpy as np
cimport numpy as np


cdef class PotentialOperator:
    cdef DiscreteBoundaryOperatorBase _op
    cdef int _component_count
    cdef Space _space
    cdef np.ndarray _evaluation_points

