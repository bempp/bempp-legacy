from bempp.core.utils cimport Matrix
from bempp.core.utils cimport eigen_matrix_to_np_float64
from libcpp.vector cimport vector


cdef extern from "bempp/fiber/numerical_quadrature.hpp" namespace "Fiber":
    cdef void c_fill_single_quadrature_points_and_weights "Fiber::fillSingleQuadraturePointsAndWeights<double>"(
              int elementCornerCount,
              int accuracyOrder,
              Matrix[double]& points,
              vector[double]& weights)


def fill_quadrature_points_and_weights(int element_corner_count, int accuracy_order):
    """Return the quadrature points and weights for the given element type and accuracy."""

    import numpy as np
    
    cdef Matrix[double] points
    cdef vector[double] weights

    c_fill_single_quadrature_points_and_weights(element_corner_count, accuracy_order,
            points, weights)

    return (eigen_matrix_to_np_float64(points), np.array(weights))



