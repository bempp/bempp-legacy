
def lu_solve(A, b):

    from numpy.linalg import solve
    from bempp import GridFunction

    mat = A.weak_form().as_matrix()
    b = b.coefficients

    sol = solve(mat,b)
    return GridFunction(A.domain,coefficients=sol) 


