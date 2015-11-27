"""Definition of a Bem++ space object and the associated factory function."""

class Space(object):
    """ Space of functions defined on a grid

        Attributes
        ----------
        grid : bempp.api.grid.Grid
            Grid over which to discretize the space.

        dtype : numpy.dtype
            Type of the basis functions in this space.

        codomain_dimension : int
            Number of components of values of functions in this space.

        domain_dimension : int
            Dimension of the domain on which the space is defined.

        global_dof_count : int
            Number of global degrees of freedom.

        flat_local_dof_count : int
            Total number of local degrees of freedom.

        global_dof_interpolation_points : np.ndarray 
            (3xN) matrix of global interpolation points for the space,
            where each column is the coordinate of an interpolation point.

        global_dof_interpolation_points : np.ndarray
            (3xN) matrix of normal directions associated with the interpolation points.

    """
    
    def __init__(self, impl):
        self._impl = impl

    def __eq__(self, other):
        return self.is_identical(other)

    def __ne__(self, other):
        return not self.is_identical(other)

    def is_compatible(self, other):
        """Return true if spaces have the same number of global dofs."""
        return self._impl.is_compatible(other._impl)

    def is_identical(self, other):
        """Return true of spaces are identical."""
        return self._impl.is_identical(other._impl)

    def get_global_dofs(self, element, dof_weights=False):
        """Return the global dofs associated with the local dofs on the given element.

        If `dof_weights=True` also return the associated dof_weights

        """
        return self._impl.get_global_dofs(element._impl, dof_weights)

    def shapeset(self, element):
        """Return the Shapeset associated with a given element."""

        from bempp.api.space.shapeset import Shapeset
        return Shapeset(
                self._impl.shapeset(element._impl))



    def evaluate_local_basis(self, element, local_coordinates,
                             local_coefficients):
        """Evaluate a local basis on a given element."""
        return self._impl.evaluate_local_basis(element._impl, 
                                               local_coordinates,
                                               local_coefficients)

    @property
    def dtype(self):
        """Return the data type of the basis functions in the space."""
        return self._impl.dtype

    @property
    def codomain_dimension(self):
        """Return the number of components of values of functions in this space (e.g. 1 for scalar functions)."""
        return self._impl.codomain_dimension

    @property
    def grid(self):
        """Return the underlying grid of the space."""
        from bempp.api.grid import Grid
        return Grid(self._impl.grid)

    @property
    def domain_dimension(self):
        """Dimension of the domain on which the space is defined."""
        return self._impl.domain_dimension

    @property
    def global_dof_count(self):
        """Return the number of global degrees of freedom."""
        return self._impl.global_dof_count

    @property
    def flat_local_dof_clount(self):
        """Return the total number of local degrees of freedom."""
        return self._impl.flat_local_dof_count

    @property
    def discontinuous_space(self):
        """Return the associated discontinuous scalar space."""
        return Space(self._impl.discontinuous_space)

    @property
    def global_dof_interpolation_points(self):
        """ Return a (3xN) matrix of the N global interpolation points for the space.
        """
        return self._impl.global_dof_interpolation_points

    @property
    def global_dof_normals(self):
        """ Return a (3xN) matrix of normal directions associated with the interpolation points.
        """
        return self._impl.global_dof_normals


def function_space(grid, kind, order, domains=None, closed=True):
    """ Return a space defined over a given grid.

    Parameters
    ----------
    grid : bempp.api.grid.Grid
        The grid object over which the space is defined.

    kind : string
        The type of space. Currently, the following types
        are supported:
        "P" : Continuous and piecewise polynomial functions.
        "DP" : Discontinuous and elementwise polynomial functions.
        "B-P": Polynomial spaces on barycentric grids.
        "B-DP": Polynomial discontinuous spaces on barycentric grids.
        "DUAL": Dual space on dual grid (only implemented for constants).
        "RT": Raviart-Thomas Vector spaces.

    order : int
        The order of the space, e.g. 0 for piecewise const, 1 for
        piecewise linear functions.

    domains : list
        List of integers specifying a list of physical entities
        of subdomains that should be included in the space.

    closed : bool
        Specifies whether the space is defined on a closed
        or open subspace.

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
    from bempp.core.space.space import function_space as _function_space
    return Space(_function_space(grid._impl, kind, order, domains, closed))


