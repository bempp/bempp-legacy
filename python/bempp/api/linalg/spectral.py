"""A collection of spectral tools useful for surface problems."""



def laplace_beltrami_eigenpairs(space, k=6, tol=1E-5):
    """Compute k Laplace Beltrami eigenpairs over a given space."""

    import bempp.api
    
    from scipy.sparse.linalg import eigs

    op = bempp.api.operators.boundary.sparse.laplace_beltrami(space, space, space)
    ident = bempp.api.operators.boundary.sparse.identity(space, space, space)

    return eigs(op, M=mass, sigma=0, tol=tol, k=k)



