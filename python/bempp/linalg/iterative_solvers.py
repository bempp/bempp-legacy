import scipy.sparse.linalg
from bempp.assembly import BoundaryOperatorBase, GridFunction

def gmres(A, b, tol=1E-5, restart=None, maxiter=None, M=None, callback=None):

    if not isinstance(A,BoundaryOperatorBase):
        raise ValueError("A must be of type BoundaryOperatorBase")

    if not isinstance(b,GridFunction):
        raise ValueError("b must be of type GridFunction")

    x, info = scipy.sparse.linalg.gmres(A.weak_form(), b.projections(A.dual_to_range), 
            tol=tol, restart=restart, maxiter=maxiter, M=M, callback=callback)

    return (GridFunction(A.domain, result_type=b.result_type,coefficients = x.ravel()),
            info)



def cg(A, b, tol=1E-5, maxiter=None, M=None, callback=None):

    if not isinstance(A,BoundaryOperatorBase):
        raise ValueError("A must be of type BoundaryOperatorBase")

    if not isinstance(b,GridFunction):
        raise ValueError("b must be of type GridFunction")

    x, info = scipy.sparse.linalg.cg(A.weak_form(), b.projections(A.dual_to_range), 
            tol=tol, maxiter=maxiter, M=M, callback=callback)

    return (GridFunction(A.domain, result_type=b.result_type,coefficients = x.ravel()),
            info)
    


