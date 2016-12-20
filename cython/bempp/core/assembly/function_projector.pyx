cimport numpy
import numpy
numpy.import_array()

from bempp.core.utils cimport Vector
from bempp.core.utils cimport complex_double
from bempp.core.utils cimport c_ParameterList
from bempp.core.utils cimport ParameterList
from bempp.core.utils cimport catch_exception
from bempp.core.space cimport Space, c_Space
from cython.operator cimport dereference as deref

cdef extern from "bempp/core/assembly/function_projector.hpp" namespace "Bempp":
    cdef object calculateProjections "Bempp::calculateProjections<double, std::complex<double>>"(
            const c_ParameterList&, object, const c_Space[double]&) except +catch_exception


def calculate_projection(ParameterList parameters not None, object fun, Space space not None):
    """Compute the projection of a Python function onto a function space."""

    import numpy as np

    res =  calculateProjections(deref(parameters.impl_), fun, deref(space.impl_))
    # Check if function projection is real. If yes, return only real part of array.
    if (np.isreal(res).all()):
        return np.real(res)
    else:
        return res






