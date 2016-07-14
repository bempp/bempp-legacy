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

    def global_to_local_dofs(self, global_dofs):
        """Return the local dofs and weights for the given list of global dofs."""

        import numpy as np
        if np.min(global_dofs) < 0 or np.max(global_dofs) >= self.global_dof_count:
            raise ValueError(
                "For each dof index i it must hold that 0 <=i < space.global_dof_count")

        return self._impl.global_to_local_dofs(global_dofs)

    def evaluate_surface_gradient(self, element, local_coordinates, local_coefficients):
        """Evaluate the local surface gradient on a given element."""

        if self.codomain_dimension > 1:
            raise ValueError("Method only implemented for scalar spaces.")

        from bempp.core.space.space import evaluate_local_surface_gradient_ext
        return evaluate_local_surface_gradient_ext(self._impl, element._impl,
                                                   local_coordinates,
                                                   local_coefficients)

    @property
    def evaluation_functor(self):
        """Returns the functor transformation of shapesets from the reference element."""

        return self._evaluation_functor

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
    def flat_local_dof_count(self):
        """Return the total number of local degrees of freedom."""
        return self._impl.flat_local_dof_count

    @property
    def discontinuous_space(self):
        """Return the associated discontinuous scalar space."""
        return self._discontinuous_space

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

    @property
    def has_non_barycentric_space(self):
        """ Return if the space has an equivalent space that lives on the original grid."""
        return self._has_non_barycentric_space

    @property
    def non_barycentric_space(self):
        return self._non_barycentric_space

    @property
    def order(self):
        return self._order


class DiscontinuousPolynomialSpace(Space):
    """Represents a space of discontinuous, polynomial functions."""

    def __init__(self, grid, order, domains=None, closed=True, reference_point_on_segment=True,
                 element_on_segment=False):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import scalar_function_value_functor

        super(DiscontinuousPolynomialSpace, self).__init__(
            _function_space(grid._impl, "DP", order, domains, closed,
                            False, reference_point_on_segment, element_on_segment))

        self._order = order
        self._has_non_barycentric_space = True
        self._non_barycentric_space = self
        self._discontinuous_space = self
        self._evaluation_functor = scalar_function_value_functor()


class BarycentricDiscontinuousPolynomialSpace(Space):
    """Represents a space of discontinuous, polynomial functions over a barycentric grid."""

    def __init__(self, grid, order):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import scalar_function_value_functor

        super(BarycentricDiscontinuousPolynomialSpace, self).__init__(
            _function_space(grid._impl, "B-DP", order))

        self._order = order
        self._has_non_barycentric_space = True
        self._non_barycentric_space = DiscontinuousPolynomialSpace(grid, order)
        self._discontinuous_space = DiscontinuousPolynomialSpace(
            grid.barycentric_grid(), order)
        self._evaluation_functor = scalar_function_value_functor()


class ContinuousPolynomialSpace(Space):
    """Represents a space of continuous, polynomial functions."""

    def __init__(self, grid, order, domains=None, closed=True, strictly_on_segment=False,
                 element_on_segment=False):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import scalar_function_value_functor

        super(ContinuousPolynomialSpace, self).__init__(
            _function_space(grid._impl, "P", order, domains, closed,
                            strictly_on_segment, True, element_on_segment))

        self._order = order
        self._has_non_barycentric_space = True
        self._non_barycentric_space = self
        self._discontinuous_space = DiscontinuousPolynomialSpace(grid, order)
        self._evaluation_functor = scalar_function_value_functor()

class BarycentricContinuousPolynomialSpace(Space):
    """Represents a space of continuous, polynomial functions on a barycentric grid."""

    def __init__(self, grid, order):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import scalar_function_value_functor

        super(BarycentricContinuousPolynomialSpace, self).__init__(
            _function_space(grid._impl, "B-P", order))

        self._order = order
        self._has_non_barycentric_space = True
        self._non_barycentric_space = ContinuousPolynomialSpace(grid, order)
        self._discontinuous_space = DiscontinuousPolynomialSpace(
            grid.barycentric_grid(), order)
        self._evaluation_functor = scalar_function_value_functor()


class DualSpace(Space):
    """A space of piecewise constant dual functions over a barycentric grid."""

    def __init__(self, grid):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import scalar_function_value_functor

        super(DualSpace, self).__init__(
            _function_space(grid._impl, "DUAL", 0))

        self._order = 0
        self._has_non_barycentric_space = False
        self._non_barycentric_space = None
        self._discontinuous_polynomial_space = DiscontinuousPolynomialSpace(
            grid.barycentric_grid(), 1)
        self._evaluation_functor = scalar_function_value_functor()


class RTSpace(Space):
    """A space of Raviart-Thomas functions."""

    def __init__(self, grid, domains, closed):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hdiv_function_value_functor

        super(RTSpace, self).__init__(
            _function_space(grid._impl, "RT", 0, domains, closed))

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = self
        self._discontinuous_space = DiscontinuousPolynomialSpace(grid, 1)
        self._evaluation_functor = hdiv_function_value_functor()

class NCSpace(Space):
    """A space of Nedelec functions."""

    def __init__(self, grid, domains, closed):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hcurl_function_value_functor

        super(NCSpace, self).__init__(
            _function_space(grid._impl, "NC", 0, domains, closed))

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = self
        self._discontinuous_space = DiscontinuousPolynomialSpace(grid, 1)
        self._evaluation_functor = hcurl_function_value_functor()
        self._hdiv_space = RTSpace(grid, domains, closed)

class RWGSpace(Space):
    """A space of RWG functions."""

    def __init__(self, grid, domains, closed):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hdiv_function_value_functor

        super(RWGSpace, self).__init__(
            _function_space(grid._impl, "RWG", 0, domains, closed))

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = self
        self._discontinuous_space = DiscontinuousPolynomialSpace(grid, 1)
        self._evaluation_functor = hdiv_function_value_functor()

class SNCSpace(Space):
    """A space of scaled Nedelec functions."""

    def __init__(self, grid, domains, closed):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hcurl_function_value_functor

        super(SNCSpace, self).__init__(
            _function_space(grid._impl, "SNC", 0, domains, closed))

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = self
        self._discontinuous_space = DiscontinuousPolynomialSpace(grid, 1)
        self._evaluation_functor = hcurl_function_value_functor()
        self._hdiv_space = RWGSpace(grid, domains, closed)

class BarycentricRTSpace(Space):
    """A space of Raviart-Thomas functions on a barycentric grid."""

    def __init__(self, grid):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hdiv_function_value_functor

        super(BarycentricRTSpace, self).__init__(
            _function_space(grid._impl, "B-RT", 0))

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = RTSpace(grid, None, True)
        self._discontinuous_space = DiscontinuousPolynomialSpace(
            grid.barycentric_grid(), 1)
        self._evaluation_functor = hdiv_function_value_functor()

class BarycentricNCSpace(Space):
    """A space of Nedelec functions on a barycentric grid."""

    def __init__(self, grid):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hcurl_function_value_functor

        super(BarycentricNCSpace, self).__init__(
            _function_space(grid._impl, "B-NC", 0))

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = NCSpace(grid, None, True)
        self._discontinuous_space = DiscontinuousPolynomialSpace(
            grid.barycentric_grid(), 1)
        self._evaluation_functor = hcurl_function_value_functor()
        self._hdiv_space = BarycentricRTSpace(grid)



class BarycentricRWGSpace(Space):
    """A space of RWG functions on a barycentric grid."""

    def __init__(self, grid):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hdiv_function_value_functor

        super(BarycentricRWGSpace, self).__init__(
            _function_space(grid._impl, "B-RWG", 0))

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = RWGSpace(grid, None, True)
        self._discontinuous_space = DiscontinuousPolynomialSpace(
            grid.barycentric_grid(), 1)
        self._evaluation_functor = hdiv_function_value_functor()

class BarycentricSNCSpace(Space):
    """A space of scaled Nedelec functions on a barycentric grid."""

    def __init__(self, grid):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hcurl_function_value_functor

        super(BarycentricSNCSpace, self).__init__(
            _function_space(grid._impl, "B-SNC", 0))

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = NCSpace(grid, None, True)
        self._discontinuous_space = DiscontinuousPolynomialSpace(
            grid.barycentric_grid(), 1)
        self._evaluation_functor = hcurl_function_value_functor()
        self._hdiv_space = BarycentricRWGSpace(grid)




class BuffaChristiansenSpace(Space):
    """A space of Buffa-Christiansen basis functions on a barycentrid grid."""

    def __init__(self, grid):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hdiv_function_value_functor

        super(BuffaChristiansenSpace, self).__init__(_function_space(
            grid._impl, "BC", 0))

        self._order = 0
        self._has_non_barycentric_space = False
        self._non_barycentric_space = None
        self._discontinuous_space = DiscontinuousPolynomialSpace(
            grid.barycentric_grid(), 1)
        self._evaluation_functor = hdiv_function_value_functor()

class RotatedBuffaChristiansenSpace(Space):
    """A space of rotated Buffa-Christiansen curl conforming basis functions on a barycentrid grid."""

    def __init__(self, grid):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hcurl_function_value_functor

        super(RotatedBuffaChristiansenSpace, self).__init__(_function_space(
            grid._impl, "RBC", 0))

        self._order = 0
        self._has_non_barycentric_space = False
        self._non_barycentric_space = None
        self._discontinuous_space = DiscontinuousPolynomialSpace(
            grid.barycentric_grid(), 1)
        self._evaluation_functor = hcurl_function_value_functor()
        self._hdiv_space = BuffaChristiansenSpace(grid)

def function_space(grid, kind, order, domains=None, closed=True, strictly_on_segment=False,
                   reference_point_on_segment=True, element_on_segment=False):
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
            "RT": Raviart-Thomas Vector spaces.
            "RWG": RWG Vector spaces.
            "NC": Nedelec Vector spaces.
            "SNC": Scaled Nedelec Vector spaces. The Nedelec basis functions
                   are scaled with the edge length so that they are identical
                   to RWG functions crossed with the element normals.

            "B-P": Polynomial spaces on barycentric grids.
            "B-DP": Polynomial discontinuous spaces on barycentric grids.
            "B-RT": Raviart-Thomas Vector spaces on barycentric grids.
            "B-RWG": RWG Vector spaces on barycentric grids.
            "B-NC": Nedelec Vector spaces on barycentric grids.
            "B-SNC": Scaled Nedelec Vector spaces on barycentric grids.

            "DUAL": Dual space on dual grid (only implemented for constants).
            "BC": Buffa-Christian Vector space.
            "RBC": Rotated Buffa-Christian Vector space of curl-conforming
                   functions.

    order : int
        The order of the space, e.g. 0 for piecewise const, 1 for
        piecewise linear functions.

    domains : list
        List of integers specifying a list of physical entities
        of subdomains that should be included in the space.

    closed : bool
        Specifies whether the space is defined on a closed
        or open subspace.

    strictly_on_segment: bool
        Specifies whether local basis functions are truncated to
        the domains specified (True) or if they are allowed to extend
        past the domains (False). Default is False. This argument is
        only used for scalar continuous spaces.

    reference_point_on_segment: bool
        If true only include a dof if its reference point (i.e. the
        dof position) is part of the segment. This argument is only
        used for discontinuous spaces (default is True).

    element_on_segment: bool
        If true restrict the dofs to those whose support element
        is part of the segment (default is False).




    Notes
    -----
    The most frequent used types are the space of piecewise constant
    functions (kind="DP", order=0) and the space of continuous,
    piecewise linear functions (kind="P", order=1).

    Either one of `reference_point_on_segment` or `element_on_segment`
    must be true for discontinuous spaces. For piecewise constant spaces
    neither of these two options has any effect.

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
    if kind == "DP":
        return DiscontinuousPolynomialSpace(grid, order, domains, closed, reference_point_on_segment,
                                            element_on_segment)
    elif kind == "B-DP":
        return BarycentricDiscontinuousPolynomialSpace(grid, order)
    elif kind == "P":
        return ContinuousPolynomialSpace(grid, order, domains, closed, strictly_on_segment,
                                         element_on_segment)
    elif kind == "B-P":
        return BarycentricContinuousPolynomialSpace(grid, order)
    elif kind == "DUAL":
        if order > 0:
            raise ValueError("Only order zero dual spaces supported.")
        return DualSpace(grid)
    elif kind == "RT":
        if order > 0:
            raise ValueError(
                "Only order zero Raviart-Thomas spaces supported.")
        return RTSpace(grid, domains, closed)
    elif kind == "NC":
        if order > 0:
            raise ValueError(
                "Only order zero Nedelec spaces supported.")
        return NCSpace(grid, domains, closed)
    elif kind == "SNC":
        if order > 0:
            raise ValueError(
                "Only order zero Nedelec spaces supported.")
        return SNCSpace(grid, domains, closed)
    elif kind == "RWG":
        if order > 0:
            raise ValueError("Only order zero RWG spaces supported.")
        return RWGSpace(grid, domains, closed)
    elif kind == "B-RT":
        if order > 0:
            raise ValueError(
                "Only order zero Raviart-Thomas functions supported.")
        return BarycentricRTSpace(grid)
    elif kind == "B-NC":
        if order > 0:
            raise ValueError(
                "Only order zero Nedelec functions supported.")
        return BarycentricNCSpace(grid)
    elif kind == "B-SNC":
        if order > 0:
            raise ValueError(
                "Only order zero Nedelec functions supported.")
        return BarycentricSNCSpace(grid)
    elif kind == "B-RWG":
        if order > 0:
            raise ValueError("Only order zero RWG functions supported.")
        return BarycentricRWGSpace(grid)
    elif kind == "BC":
        if order > 0:
            raise ValueError(
                "Only order zero Buffa-Christiansen functions supported.")
        return BuffaChristiansenSpace(grid)
    elif kind == "RBC":
        if order > 0:
            raise ValueError(
                "Only order zero Buffa-Christiansen functions supported.")
        return RotatedBuffaChristiansenSpace(grid)
    else:
        raise ValueError("Unknown space type.")

