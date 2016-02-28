import scipy.sparse.linalg
import numpy as np
from bempp.api.assembly import GridFunction
from bempp.api.assembly import BoundaryOperator

class _it_counter(object):

    def __init__(self, store_residuals):
        self._count = 0
        self._store_residuals = store_residuals
        self._residuals = []

    def __call__(self, x):
        self._count += 1
        if self._store_residuals:
            self._residuals.append(np.linalg.norm(x))


    @property
    def count(self):
        return self._count

    @property
    def residuals(self):
        return self._residuals

def gmres(A, b, tol=1E-5, restart=None, maxiter=None, use_strong_form=False, return_residuals=False):
    """Interface to the scipy.sparse.linalg.gmres function.

    This function behaves like the scipy.sparse.linalg.gmres function. But
    instead of a linear operator and a vector b it takes a boundary operator
    and a grid function. The result is returned as a grid function in the
    correct space.

    """
    import bempp.api
    import time

    if not isinstance(A, BoundaryOperator):
        raise ValueError("A must be of type BoundaryOperator")

    if not isinstance(b, GridFunction):
        raise ValueError("b must be of type GridFunction")

    # Assemble weak form before the logging messages

    if use_strong_form:
        if not A.range.is_compatible(b.space):
            raise ValueError("The range of A and the space of A must have the same number of unknowns if the strong form is used.")
        A_op = A.strong_form()
        b_vec = b.coefficients
    else:
        A_op = A.weak_form()
        b_vec = b.projections(A.dual_to_range)

    callback = _it_counter(return_residuals)

    bempp.api.LOGGER.info("Starting GMRES iteration")
    start_time = time.time()
    x, info = scipy.sparse.linalg.gmres(A_op, b_vec,
                                        tol=tol, restart=restart, maxiter=maxiter, callback=callback)
    end_time = time.time()
    bempp.api.LOGGER.info("GMRES finished in {0} iterations and took {1:.2E} sec.".format(
        callback.count, end_time - start_time))

    res_fun = GridFunction(A.domain, coefficients=x.ravel())

    if return_residuals:
        return res_fun, info, callback.residuals
    else:
        return res_fun, info


def cg(A, b, tol=1E-5, maxiter=None, 
        use_strong_form=False, return_residuals=False):
    """Interface to the scipy.sparse.linalg.cg function.

    This function behaves like the scipy.sparse.linalg.cg function. But
    instead of a linear operator and a vector b it takes a boundary operator
    and a grid function. The result is returned as a grid function in the
    correct space.

    """
    import bempp.api
    import time

    if not isinstance(A, BoundaryOperator):
        raise ValueError("A must be of type BoundaryOperator")

    if not isinstance(b, GridFunction):
        raise ValueError("b must be of type GridFunction")

    if use_strong_form:
        if not A.range.is_compatible(b.space):
            raise ValueError("The range of A and the space of A must have the same number of unknowns if the strong form is used.")
        A_op = A.strong_form()
        b_vec = b.coefficients
    else:
        A_op = A.weak_form()
        b_vec = b.projections(A.dual_to_range)

    callback = _it_counter(return_residuals)
    bempp.api.LOGGER.info("Starting CG iteration")
    start_time = time.time()
    x, info = scipy.sparse.linalg.cg(A_op, b_vec,
                                     tol=tol, maxiter=maxiter, callback=callback)
    end_time = time.time()
    bempp.api.LOGGER.info("CG finished in {0} iterations and took {1:.2E} sec.".format(
        callback.count, end_time - start_time))

    res_fun = GridFunction(A.domain, coefficients=x.ravel())

    if return_residuals:
        return res_fun, info, callback.residuals
    else:
        return res_fun, info
