import scipy.sparse.linalg
from bempp.assembly import BoundaryOperator, GridFunction
from numpy import ndarray as _ndarray

def gmres(A, b, tol=1E-5, restart=None, maxiter=None, M=None, callback=None):

    if not isinstance(A,BoundaryOperator):
        raise ValueError("A must be of type BoundaryOperator")

    if isinstance(b,GridFunction):
        b_arr = b.projections(A.dual_to_range)
    elif isinstance(b,_ndarray):
        b_arr = b
    else:
        raise ValueError("b must be of type GridFunction or numpy.array")

    x, info = scipy.sparse.linalg.gmres(A.weak_form(), b_arr,
            tol=tol, restart=restart, maxiter=maxiter, M=M, callback=callback)

    return (GridFunction(A.domain, coefficients=x.ravel()), info)



def cg(A, b, tol=1E-5, maxiter=None, M=None, callback=None):

    if not isinstance(A,BoundaryOperator):
        raise ValueError("A must be of type BoundaryOperator")

    if isinstance(b,GridFunction):
        b_arr = b.projections(A.dual_to_range)
    elif isinstance(b,_ndarray):
        b_arr = b
    else:
        raise ValueError("b must be of type GridFunction or numpy.array")

    x, info = scipy.sparse.linalg.cg(A.weak_form(), b_arr,
            tol=tol, maxiter=maxiter, M=M, callback=callback)

    return (GridFunction(A.domain, coefficients=x.ravel()), info)
    


