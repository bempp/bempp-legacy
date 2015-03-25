
def lu(A, b):

    from numpy.linalg import solve
    from bempp import GridFunction

    mat = A.weak_form().as_matrix()
    vec = b.projections(A.dual_to_range)

    sol = solve(mat,vec)
    return GridFunction(A.domain,coefficients=sol) 


