from .space cimport c_Space, _py_get_space_ptr, SpaceVariants
from . import space as s


def space(grid, kind, order, dtype='float64'):

    if kind=="P":
        if not (order>=1 and order <=10):
            raise ValueError("Order must be between 1 and 10")
        return s.PiecewisePolynomialContinuousScalarSpace(grid,dtype,order)
    elif kind=="DP":
        return s.PiecewisePolynomialDiscontinuousScalarSpace(grid,dtype,order)
    else:
        raise ValueError("Unknown kind")

