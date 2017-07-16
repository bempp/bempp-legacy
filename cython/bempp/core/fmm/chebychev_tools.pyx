from cython.operator cimport dereference as deref
from bempp.core.utils.shared_ptr cimport shared_ptr
from bempp.core.utils.eigen cimport Vector
from bempp.core.utils.eigen cimport eigen_vector_to_np_float64
from bempp.core.utils.eigen cimport np_to_eigen_vector_float64
from bempp.core.utils.eigen cimport eigen_matrix_to_np_float64


cdef class ChebychevTools:

    def __cinit__(self, int order):
        pass

    def __init__(self, int order):
        if order == 0:
            raise ValueError("'order' must at least be 1")
        self.impl_.assign(shared_ptr[c_ChebychevTools](new c_ChebychevTools(order)))

    def __dealloc__(self):
        self.impl_.reset()

    def chebychev_nodes(self):
        """Compute the Chebychev nodes"""
        return eigen_vector_to_np_float64(deref(self.impl_).chebychevNodes())

    def chebychev_values(self):
        """Return the values of the Cheb. Pol. at the nodes."""
        return eigen_matrix_to_np_float64(deref(self.impl_).chebychevPolValuesAtNodes())

    def evaluate_interpolation_polynomial(self, weights, evaluation_points):
        """Evaluate an interp. polynomial with given weights at the given points."""

        cdef Vector[double] result
        deref(self.impl_).evaluateInterpolationPolynomial(
            np_to_eigen_vector_float64(weights),
            np_to_eigen_vector_float64(evaluation_points),
            result)
        return eigen_vector_to_np_float64(result)

    def interpolate_to_children(self, parent_length, child_length):
        """Return inteprolation matrix to children"""
        return eigen_matrix_to_np_float64(
                deref(self.impl_).interpolateToChildren(parent_length, child_length))







