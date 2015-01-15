__all__=['function_space']
from . import space

def function_space(grid, kind, order, dtype='float64'):
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

    dtype : np.dtype, optional
        The data type of values of the basis functions in the space.
        The default is 'np.float64'.

    Notes
    -----
    The most frequent used types are the space of piecewise constant
    functions (kind="DP", order=0) and the space of continuous,
    piecewise linear functions (kind="DP", order=0).

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
            s = space.PiecewiseLinearContinuousScalarSpace(grid,order,dtype)
        else:
            s = space.PiecewisePolynomialContinuousScalarSpace(grid,order,dtype)
    elif kind=="DP":
        if not (order>=0 and order <=10):
            raise ValueError("Order must be between 0 and 10")
        if (order==0):
            s = space.PiecewiseConstantScalarSpace(grid,order,dtype)
        else:
            s = space.PiecewisePolynomialDiscontinuousScalarSpace(grid,order,dtype)
    else:
        raise ValueError("Unknown kind")

    return s

