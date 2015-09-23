import scipy.sparse.linalg
from bempp.api.assembly import GridFunction
from bempp.api.assembly import BoundaryOperator


def gmres(A, b, tol=1E-5, restart=None, maxiter=None, M=None, callback=None):
    if not isinstance(A, BoundaryOperator):
        raise ValueError("A must be of type BoundaryOperator")

    if not isinstance(b, GridFunction):
        raise ValueError("b must be of type GridFunction")

    x, info = scipy.sparse.linalg.gmres(A.weak_form(), b.projections(A.dual_to_range),
                                        tol=tol, restart=restart, maxiter=maxiter, M=M, callback=callback)

    return GridFunction(A.range, coefficients=x.ravel()), info


def cg(A, b, tol=1E-5, maxiter=None, M=None, callback=None):
    if not isinstance(A, BoundaryOperator):
        raise ValueError("A must be of type BoundaryOperator")

    if not isinstance(b, GridFunction):
        raise ValueError("b must be of type GridFunction")

    x, info = scipy.sparse.linalg.cg(A.weak_form(), b.projections(A.dual_to_range),
                                     tol=tol, maxiter=maxiter, M=M, callback=callback)

    return GridFunction(A.range, coefficients=x.ravel()), info
