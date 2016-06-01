
def lu(A, b):
    """Simple direct solver interface.

    This function takes an operator and a grid function,
    converts the operator into a dense matrix and solves
    the system via LU decomposition. The result is again
    returned as a grid function.

    """

    from numpy.linalg import solve
    from bempp.api import GridFunction, as_matrix

    mat = as_matrix(A.weak_form())
    vec = b.projections(A.dual_to_range)

    sol = solve(mat, vec)
    return GridFunction(A.domain, coefficients=sol)
