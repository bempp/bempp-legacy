<%
from data_types import dtypes
%>

from bempp.utils.armadillo cimport Mat
from bempp.utils.enum_types cimport transposition_mode
from bempp.utils cimport complex_float,complex_double

cdef class DiscreteBoundaryOperator:

    def __cinit__(self):
        pass

