

def gauss_triangle_points_and_weights(accuracy_order):
    """Return Gauss integration rules over triangles.

    This function returns the points and weights for
    the symmetric Gauss rule over the reference triangle.

    """

    from bempp.core.fiber.integrators import \
            fill_quadrature_points_and_weights

    return fill_quadrature_points_and_weights(3, accuracy_order)
