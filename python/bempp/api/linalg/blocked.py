import scipy.sparse.linalg
from bempp.api.utils.blocked_linear_operator import BlockedDiscreteLinearOperator

def gmres(A, b, tol=1E-5, restart=None, maxiter=None, M=None, callback=None):

    if not isinstance(A,BlockedDiscreteLinearOperator):
        raise ValueError("A must be of type BlockedDiscreteLinearOperator")

    soln,info = scipy.sparse.linalg.gmres(A, b, tol=tol, restart=restart, maxiter=maxiter, M=M, callback=callback)
    output = []
    curr = 0
    for column_size in A.column_dimensions:
        prev = curr
        curr += column_size
        output.append(soln[prev:curr])

    return output,info


def cg(A, b, tol=1E-5, maxiter=None, M=None, callback=None):

    if not isinstance(A,BlockedDiscreteLinearOperator):
        raise ValueError("A must be of type BlockedDiscreteLinearOperator")

    soln,info = scipy.sparse.linalg.cg(A, b, tol=tol, restart=restart, maxiter=maxiter, M=M, callback=callback)
    output = []
    curr = 0
    for column_size in A.column_dimensions:
        prev = curr
        curr += column_size
        output.append(soln[prev:curr])

    return output,info
