__all__=['function_space']
from . import space

def function_space(grid, kind, order):
    """ 

    Return a space defined over a given grid.

    Parameters
    ----------
    grid : bempp.Grid
        The grid object over which the space is defined.

    kind : string
        The type of space. Currently, the following types
        are supported:
        "P" : Continuous and piecewise polynomial functions.
        "DP" : Discontinuous and elementwise polynomial functions.

    order : int
        The order of the space, e.g. 0 for piecewise const, 1 for
        piecewise linear functions.

    Notes
    -----
    The most frequent used types are the space of piecewise constant
    functions (kind="DP", order=0) and the space of continuous,
    piecewise linear functions (kind="P", order=1).

    This is a factory function that initializes a space object. To 
    see a detailed help for space objects see the documentation
    of the instantiated object.

    Examples
    --------
    To initialize a space of piecewise constant functions use

    >>> space = function_space(grid,"DP",0)

    To initialize a space of continuous, piecewise linear functions, use

    >>> space = function_space(grid,"P",1)

    """


    if kind=="P":
        if not (order>=1 and order <=10):
            raise ValueError("Order must be between 1 and 10")
        if (order==1):
            s = space.PiecewiseLinearContinuousScalarSpace(grid,order)
        else:
            s = space.PiecewisePolynomialContinuousScalarSpace(grid,order)
    elif kind=="DP":
        if not (order>=0 and order <=10):
            raise ValueError("Order must be between 0 and 10")
        if (order==0):
            s = space.PiecewiseConstantScalarSpace(grid,order)
        else:
            s = space.PiecewisePolynomialDiscontinuousScalarSpace(grid,order)

    elif kind=="RT":
        if order!=0:
            raise ValueError("Only 0 order Raviart-Thomas spaces are implemented.")
        s = space.RaviartThomas0VectorSpace(grid,order)
    else:
        raise ValueError("Unknown kind")

    return s

