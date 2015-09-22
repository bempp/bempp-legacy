
def lu(A, b):

    from numpy.linalg import solve
    from bempp.api import GridFunction, as_matrix

    mat = as_matrix(A.weak_form())
    vec = b.projections(A.dual_to_range)

    sol = solve(mat,vec)
    return GridFunction(A.domain,coefficients=sol) 


