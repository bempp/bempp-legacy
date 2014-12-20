from bempp.utils.shared_ptr cimport shared_ptr
from bempp.utils.unique_ptr cimport unique_ptr
from bempp.utils.parameter_list cimport ParameterList

cdef extern from "bempp/utils/py_utils.hpp" namespace "Bempp":
    void catch_exception()

# Declares complex type explicitly.
# Cython 0.20 will fail if templates are nested more than three-deep,
# as in shared_ptr[ c_Space[ complex[float] ] ]
cdef extern from "bempp/utils/py_types.hpp":
    ctypedef struct complex_float
    ctypedef struct complex_double

