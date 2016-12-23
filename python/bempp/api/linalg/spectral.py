"""A collection of spectral tools useful for surface problems."""


def laplace_beltrami_eigenpairs(space, k=6, tol=1E-5):
    """Compute k Laplace Beltrami eigenpairs over a given space."""
    from scipy.sparse.linalg import eigs
    from bempp.api.operators.boundary.sparse import \
        laplace_beltrami
    from bempp.api.operators.boundary.sparse import \
        identity

    operator = laplace_beltrami(space, space, space)
    ident = identity(space, space, space)

    return eigs(operator, M=ident, sigma=0, tol=tol, k=k)
