"""Definition of a Bem++ space object and the associated factory function."""

#pylint: disable=protected-access
#pylint: disable=no-member
#pylint: disable=too-many-instance-attributes
#pylint: disable=too-many-arguments

#pylint: disable=too-many-public-methods
#pylint: disable=no-name-in-module

class Space(object):
    """
    Space of functions defined on a grid.

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
        (3xN) matrix of normal directions associated with the interpolation
        points.

    """
    def __init__(self, impl, comp_key):
        """
        Construct a space object. Should not be called directly.

        To construct a space use the function bempp.api.function_space.

        """
        self._impl = impl
        self._global_to_local_dofs = None
        self._mass_matrix = None
        self._inverse_mass_matrix = None
        self._comp_key = comp_key

    def __eq__(self, other):
        return self._comp_key == other._comp_key

    def __ne__(self, other):
        return self._comp_key != other._comp_key

    def is_compatible(self, other):
        """Return true if spaces have the same number of global dofs."""
        return self._impl.is_compatible(other._impl)

    def is_identical(self, other):
        """Return true of spaces are identical."""
        return self == other

    def get_global_dofs(self, element, dof_weights=False):
        """Global dofs associated with the local dofs on the given element.

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
        """Local dofs and weights for the given list of global dofs."""

        if self._global_to_local_dofs is None:
            self._global_to_local_dofs = self._impl.global_to_local_dofs(
                range(self.global_dof_count))

        dofs = [self._global_to_local_dofs[0][index] for index in global_dofs]
        weights = [self._global_to_local_dofs[1][index]
                   for index in global_dofs]

        return dofs, weights

    def evaluate_surface_gradient(
            self,
            element,
            local_coordinates,
            local_coefficients):
        """Evaluate the local surface gradient on a given element."""

        if self.codomain_dimension > 1:
            raise ValueError("Method only implemented for scalar spaces.")

        #pylint: disable=no-name-in-module
        from bempp.core.space.space import evaluate_local_surface_gradient_ext
        return evaluate_local_surface_gradient_ext(self._impl, element._impl,
                                                   local_coordinates,
                                                   local_coefficients)

    def mass_matrix(self):
        """Return the mass matrix associated with this space."""

        if self._mass_matrix is None:
            from bempp.api.operators.boundary.sparse import identity
            self._mass_matrix = identity(self, self, self).weak_form()

        return self._mass_matrix

    def inverse_mass_matrix(self):
        """Return the inverse mass matrix for this space."""

        from bempp.api.assembly.discrete_boundary_operator import \
            InverseSparseDiscreteBoundaryOperator

        if self._inverse_mass_matrix is None:
            self._inverse_mass_matrix = InverseSparseDiscreteBoundaryOperator(
                self.mass_matrix())
        return self._inverse_mass_matrix

    def _get_included_faces(self, grid):
        faces = []
        for e in grid.leaf_view.entity_iterator(0):
            if len([dof for dof in self.get_global_dofs(e) if dof>=0]) > 0:
                faces.append(grid.element_insertion_index(e))
        return faces

    @property
    def evaluation_functor(self):
        """Functor transformation of shapesets from the reference element."""
        return self._evaluation_functor

    @property
    def dtype(self):
        """Return the data type of the basis functions in the space."""
        return self._impl.dtype

    @property
    def codomain_dimension(self):
        """
        Number of components of values of functions in this space.

        For a scalar space this is one. For a vector space it is three.
        """
        return self._impl.codomain_dimension

    @property
    def grid(self):
        """Return the underlying grid of the space."""
        return self._grid

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
        """(3xN) matrix of the N global interpolation points for the space."""
        return self._impl.global_dof_interpolation_points

    @property
    def global_dof_normals(self):
        """(3xN) matrix of normal directions at interp. points."""
        return self._impl.global_dof_normals

    @property
    def has_non_barycentric_space(self):
        """True if space has an equivalent non-barycentric space."""
        return self._has_non_barycentric_space

    @property
    def is_barycentric(self):
        """ True if the space is defined over a barycentric refinement."""
        return self._is_barycentric

    @property
    def non_barycentric_space(self):
        """The associated non-barycentric space if it exists."""
        return self._non_barycentric_space

    @property
    def order(self):
        """Return the order of the space."""
        return self._order

    @property
    def super_space(self):
        """
        A super space that contains the space as subspace.

        This attribute is useful for the compressed assembly of operators.
        The super space is typically a more local space (often discontinuous,
        elementwise defined) that is easier to compress than the actual space.
        The actual space is then obtained via projection.

        """
        return self._super_space

    @property
    def has_local_support(self):
        """True if support of each basis fct. is restricted to one element."""
        return self == self.discontinuous_space


class DiscontinuousPolynomialSpace(Space):
    """Represents a space of discontinuous, polynomial functions."""

    def __init__(
            self,
            grid,
            order,
            domains=None,
            closed=True,
            reference_point_on_segment=True,
            element_on_segment=False,
            comp_key=None):
        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import scalar_function_value_functor

        super(
            DiscontinuousPolynomialSpace,
            self).__init__(
                _function_space(
                    grid._impl,
                    "DP",
                    order,
                    domains,
                    closed,
                    False,
                    reference_point_on_segment,
                    element_on_segment),
                comp_key=comp_key)

        self._order = order
        self._has_non_barycentric_space = True
        self._non_barycentric_space = self
        self._discontinuous_space = self
        self._super_space = self
        self._evaluation_functor = scalar_function_value_functor()
        self._is_barycentric = False
        self._grid = grid


class CustomDPSpace(Space):
    """Represents a space of discontinuous, polynomial functions on a custom segment."""

    def __init__(
            self,
            grid,
            order,
            faces_to_include,
            nodes_to_include,
            comp_key=None):
        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import scalar_function_value_functor

        super(CustomDPSpace, self).__init__(
                _function_space(
                    grid._impl, "DP-CUSTOM", order,
                    faces_to_include=faces_to_include,
                    nodes_to_include=nodes_to_include),
                comp_key=comp_key)

        self._order = order
        self._has_non_barycentric_space = True
        self._non_barycentric_space = self
        self._discontinuous_space = self
        self._super_space = self
        self._evaluation_functor = scalar_function_value_functor()
        self._is_barycentric = False
        self._grid = grid


class BarycentricDiscontinuousPolynomialSpace(Space):
    """Represents a space of disc., polynomial fct. over a barycentric grid."""

    def __init__(self, grid, order, domains, comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import scalar_function_value_functor

        super(BarycentricDiscontinuousPolynomialSpace, self).__init__(
            _function_space(grid._impl, "B-DP", order, domains), comp_key)

        self._order = order
        self._has_non_barycentric_space = True
        self._non_barycentric_space = function_space(grid, "DP", order)
        if domains is None:
            self._discontinuous_space = function_space(grid.barycentric_grid(), "DP", order)
        else:
            faces = self._get_included_faces(grid.barycentric_grid())
            self._discontinuous_space = function_space(
                grid.barycentric_grid(), "DP-CUSTOM", order, faces_to_include=faces)
        self._super_space = self._discontinuous_space
        self._evaluation_functor = scalar_function_value_functor()
        self._is_barycentric = True
        self._grid = grid.barycentric_grid()


class ContinuousPolynomialSpace(Space):
    """Represents a space of continuous, polynomial functions."""

    def __init__(
            self,
            grid,
            order,
            domains=None,
            closed=True,
            strictly_on_segment=False,
            element_on_segment=False,
            comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import scalar_function_value_functor

        super(
            ContinuousPolynomialSpace,
            self).__init__(
                _function_space(
                    grid._impl,
                    "P",
                    order,
                    domains,
                    closed,
                    strictly_on_segment,
                    True,
                    element_on_segment),
                comp_key)

        self._order = order
        self._has_non_barycentric_space = True
        self._non_barycentric_space = self
        if not closed:
            self._discontinuous_space = function_space(
                grid,
                "DP",
                order,
                domains=domains,
                closed=closed,
                reference_point_on_segment=False,
                element_on_segment=True)
        else:
            self._discontinuous_space = function_space(
                grid,
                "DP",
                order,
                domains=domains,
                closed=closed,
                reference_point_on_segment=True,
                element_on_segment=strictly_on_segment)

        self._super_space = self._discontinuous_space
        self._evaluation_functor = scalar_function_value_functor()
        self._is_barycentric = False
        self._grid = grid


class BarycentricContinuousPolynomialSpace(Space):
    """Represents a space of cont., polynomial fct. on a barycentric grid."""

    def __init__(self, grid, order, domains, closed, strictly_on_segment, comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import scalar_function_value_functor

        super(BarycentricContinuousPolynomialSpace, self).__init__(
            _function_space(grid._impl, "B-P", order, domains, closed, strictly_on_segment), comp_key)

        self._order = order
        self._has_non_barycentric_space = True
        self._non_barycentric_space = function_space(grid, "P", order)
        if domains is None:
            self._discontinuous_space = function_space(grid.barycentric_grid(), "DP", order)
        else:
            faces = self._get_included_faces(grid.barycentric_grid())
            self._discontinuous_space = function_space(
                grid.barycentric_grid(), "DP-CUSTOM", order, faces_to_include=faces)
        self._super_space = self._discontinuous_space
        self._evaluation_functor = scalar_function_value_functor()
        self._is_barycentric = True
        self._grid = grid.barycentric_grid()


class DualConstantSpace(Space):
    """A space of piecewise constant dual functions over a barycentric grid."""

    def __init__(self, grid,
            domains=None,
            closed=True,
            strictly_on_segment=False,
            comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import scalar_function_value_functor

        super(DualConstantSpace, self).__init__(
            _function_space(grid._impl, "DUAL", 0,
                            domains=domains, closed=closed, strictly_on_segment=strictly_on_segment), comp_key)

        self._order = 0
        self._has_non_barycentric_space = False
        self._non_barycentric_space = None
        if domains is None:
            self._discontinuous_space = function_space(
                grid.barycentric_grid(), "DP", 0)
        else:
            faces = self._get_included_faces(grid.barycentric_grid())
            self._discontinuous_space = function_space(
                grid.barycentric_grid(), "DP-CUSTOM", 0, faces_to_include=faces)
        self._super_space = self._discontinuous_space
        self._evaluation_functor = scalar_function_value_functor()
        self._is_barycentric = True
        self._grid = grid.barycentric_grid()


class DualLinearSpace(Space):
    """A space of piecewise linear dual functions over a barycentric grid."""

    def __init__(self, grid,
            domains=None,
            closed=True,
            strictly_on_segment=False,
            comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import scalar_function_value_functor

        super(DualLinearSpace, self).__init__(
            _function_space(grid._impl, "DUAL", 1,
                            domains=domains, closed=closed, strictly_on_segment=strictly_on_segment), comp_key)

        self._order = 1
        self._has_non_barycentric_space = False
        self._non_barycentric_space = None
        if domains is None:
            self._discontinuous_space = function_space(
                grid.barycentric_grid(), "DP", 1)
        else:
            faces = self._get_included_faces(grid.barycentric_grid())
            self._discontinuous_space = function_space(
                grid.barycentric_grid(), "DP-CUSTOM", 1, faces_to_include=faces)
        self._super_space = self._discontinuous_space
        self._evaluation_functor = scalar_function_value_functor()
        self._is_barycentric = True
        self._grid = grid.barycentric_grid()


class RTSpace(Space):
    """A space of Raviart-Thomas functions."""

    def __init__(self, grid, domains, closed, strictly_on_segment=False, comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hdiv_function_value_functor

        super(RTSpace, self).__init__(
            _function_space(grid._impl, "RT", 0, domains, closed, strictly_on_segment=strictly_on_segment), comp_key)

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = self
        if domains is None:
            self._discontinuous_space = function_space(grid, "DP", 1)
        else:
            faces = self._get_included_faces(grid)
            self._discontinuous_space = function_space(
                grid, "DP-CUSTOM", 1, faces_to_include=faces)

        self._super_space = self
        self._evaluation_functor = hdiv_function_value_functor()
        self._is_barycentric = False
        self._grid = grid


class NCSpace(Space):
    """A space of Nedelec functions."""

    def __init__(self, grid, domains, closed, strictly_on_segment=False, comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hcurl_function_value_functor

        super(NCSpace, self).__init__(
            _function_space(grid._impl, "NC", 0, domains, closed, strictly_on_segment), comp_key)

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = self
        if domains is None:
            self._discontinuous_space = function_space(
                grid,
                "DP",
                1)
        else:
            faces = self._get_included_faces(grid)
            self._discontinuous_space = function_space(
                grid, "DP-CUSTOM", 1, faces_to_include=faces)
        self._evaluation_functor = hcurl_function_value_functor()
        self._super_space = self
        self._hdiv_space = function_space(
            grid, "RT", 0, domains=domains, closed=closed)
        self._is_barycentric = False
        self._grid = grid


class RWGSpace(Space):
    """A space of RWG functions."""

    def __init__(self, grid, domains, closed, strictly_on_segment=False, comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hdiv_function_value_functor

        super(RWGSpace, self).__init__(
            _function_space(grid._impl, "RWG", 0, domains, closed, strictly_on_segment), comp_key)

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = self
        if domains is None:
            self._discontinuous_space = function_space(grid, "DP", 1)
        else:
            faces = self._get_included_faces(grid)
            self._discontinuous_space = function_space(
                grid, "DP-CUSTOM", 1, faces_to_include=faces)
        self._super_space = self
        self._evaluation_functor = hdiv_function_value_functor()
        self._is_barycentric = False
        self._grid = grid


class SNCSpace(Space):
    """A space of scaled Nedelec functions."""

    def __init__(self, grid, domains, closed, strictly_on_segment=False, comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hcurl_function_value_functor

        super(SNCSpace, self).__init__(
            _function_space(grid._impl, "SNC", 0, domains, closed, strictly_on_segment), comp_key)

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = self
        if domains is None:
            self._discontinuous_space = function_space(grid, "DP", 1)
        else:
            faces = self._get_included_faces(grid)
            self._discontinuous_space = function_space(
                grid, "DP-CUSTOM", 1, faces_to_include=faces)
        self._evaluation_functor = hcurl_function_value_functor()
        self._super_space = self
        self._hdiv_space = function_space(grid, "RWG", 0, domains, closed)
        self._is_barycentric = False
        self._grid = grid


class BarycentricRTSpace(Space):
    """A space of Raviart-Thomas functions on a barycentric grid."""

    def __init__(self, grid, domains, closed, strictly_on_segment=False, comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hdiv_function_value_functor

        super(BarycentricRTSpace, self).__init__(
            _function_space(grid._impl, "B-RT", 0, domains, closed, strictly_on_segment=strictly_on_segment), comp_key)

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = function_space(grid, "RT", 0)
        if domains is None:
            self._discontinuous_space = function_space(grid.barycentric_grid(), "DP", 1)
        else:
            faces = self._get_included_faces(grid.barycentric_grid())
            self._discontinuous_space = function_space(
                grid.barycentric_grid(), "DP-CUSTOM", 1, faces_to_include=faces)
        self._super_space = function_space(grid.barycentric_grid(), "RT", 0)
        self._evaluation_functor = hdiv_function_value_functor()
        self._is_barycentric = True
        self._grid = grid.barycentric_grid()


class BarycentricNCSpace(Space):
    """A space of Nedelec functions on a barycentric grid."""

    def __init__(self, grid, domains, closed, strictly_on_segment=False, comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hcurl_function_value_functor

        super(BarycentricNCSpace, self).__init__(
            _function_space(grid._impl, "B-NC", 0, domains, closed, strictly_on_segment), comp_key)

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = function_space(grid, "NC", 0)
        if domains is None:
            self._discontinuous_space = function_space(grid.barycentric_grid(), "DP", 1)
        else:
            faces = self._get_included_faces(grid.barycentric_grid())
            self._discontinuous_space = function_space(
                grid.barycentric_grid(), "DP-CUSTOM", 1, faces_to_include=faces)
        self._evaluation_functor = hcurl_function_value_functor()
        self._super_space = function_space(grid.barycentric_grid(), "NC", 0)
        self._hdiv_space = function_space(grid, "B-RT", 0)
        self._is_barycentric = True
        self._grid = grid.barycentric_grid()


class BarycentricRWGSpace(Space):
    """A space of RWG functions on a barycentric grid."""

    def __init__(self, grid, domains, closed, strictly_on_segment=False, comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hdiv_function_value_functor

        super(BarycentricRWGSpace, self).__init__(
            _function_space(grid._impl, "B-RWG", 0, domains, closed, strictly_on_segment), comp_key)

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = function_space(grid, "RWG", 0)
        if domains is None:
            self._discontinuous_space = function_space(grid.barycentric_grid(), "DP", 1)
        else:
            faces = self._get_included_faces(grid.barycentric_grid())
            self._discontinuous_space = function_space(
                grid.barycentric_grid(), "DP-CUSTOM", 1, faces_to_include=faces)
        self._super_space = function_space(grid.barycentric_grid(), "RWG", 0)
        self._evaluation_functor = hdiv_function_value_functor()
        self._is_barycentric = True
        self._grid = grid.barycentric_grid()


class BarycentricSNCSpace(Space):
    """A space of scaled Nedelec functions on a barycentric grid."""

    def __init__(self, grid, domains, closed, strictly_on_segment=False, comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hcurl_function_value_functor

        super(BarycentricSNCSpace, self).__init__(
            _function_space(grid._impl, "B-SNC", 0, domains, closed, strictly_on_segment), comp_key)

        self._order = 0
        self._has_non_barycentric_space = True
        self._non_barycentric_space = function_space(grid, "SNC", 0)
        if domains is None:
            self._discontinuous_space = function_space(grid.barycentric_grid(), "DP", 1)
        else:
            faces = self._get_included_faces(grid.barycentric_grid())
            self._discontinuous_space = function_space(
                grid.barycentric_grid(), "DP-CUSTOM", 1, faces_to_include=faces)
        self._evaluation_functor = hcurl_function_value_functor()
        self._super_space = function_space(grid.barycentric_grid(), "SNC", 0)
        self._hdiv_space = function_space(grid, "B-RWG", 0)
        self._is_barycentric = True
        self._grid = grid.barycentric_grid()


class BuffaChristiansenSpace(Space):
    """A space of Buffa-Christiansen basis functions on a barycentric grid."""

    def __init__(self, grid,
                 domains=None,
                 closed=True,
                 strictly_on_segment=False,
                 comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hdiv_function_value_functor

        super(BuffaChristiansenSpace, self).__init__(
            _function_space(grid._impl, "BC", 0,
                            domains=domains, closed=closed, strictly_on_segment=strictly_on_segment), comp_key)

        self._order = 0
        self._has_non_barycentric_space = False
        self._non_barycentric_space = None
        if domains is None:
            self._discontinuous_space = function_space(grid.barycentric_grid(), "DP", 1)
        else:
            faces = self._get_included_faces(grid.barycentric_grid())
            self._discontinuous_space = function_space(
                grid.barycentric_grid(), "DP-CUSTOM", 1, faces_to_include=faces)
        self._super_space = function_space(grid.barycentric_grid(), "RWG", 0)
        self._evaluation_functor = hdiv_function_value_functor()
        self._is_barycentric = True
        self._grid = grid.barycentric_grid()


class RotatedBuffaChristiansenSpace(Space):
    """A space of rotated Buffa-Christiansen curl conforming basis functions."""

    def __init__(self, grid,
                 domains=None,
                 closed=True,
                 strictly_on_segment=False,
                 comp_key=None):

        from bempp.core.space.space import function_space as _function_space
        from bempp.api.assembly.functors import hcurl_function_value_functor

        super(RotatedBuffaChristiansenSpace, self).__init__(
            _function_space(grid._impl, "RBC", 0,
                            domains=domains, closed=closed, strictly_on_segment=strictly_on_segment), comp_key)

        self._order = 0
        self._has_non_barycentric_space = False
        self._non_barycentric_space = None
        if domains is None:
            self._discontinuous_space = function_space(grid.barycentric_grid(), "DP", 1)
        else:
            faces = self._get_included_faces(grid.barycentric_grid())
            self._discontinuous_space = function_space(
                grid.barycentric_grid(), "DP-CUSTOM", 1, faces_to_include=faces)
        self._evaluation_functor = hcurl_function_value_functor()
        self._super_space = function_space(grid.barycentric_grid(), "SNC", 0)
        self._hdiv_space = BuffaChristiansenSpace(grid)
        self._is_barycentric = True
        self._grid = grid.barycentric_grid()


#pylint: disable=too-many-return-statements
#pylint: disable=too-many-branches
def function_space(
        grid,
        kind,
        order,
        domains=None,
        closed=True,
        strictly_on_segment=False,
        reference_point_on_segment=True,
        element_on_segment=False,
        faces_to_include=None,
        nodes_to_include=None):
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

            "DP-CUSTOM" : Used internally

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

    comp_key = (
        grid,
        kind,
        order,
        None if domains is None else tuple(
            sorted(domains)),
        closed,
        strictly_on_segment,
        reference_point_on_segment,
        element_on_segment)

    if kind == "DP":
        return DiscontinuousPolynomialSpace(
            grid,
            order,
            domains,
            closed,
            reference_point_on_segment,
            element_on_segment,
            comp_key)
    elif kind == "B-DP":
        return BarycentricDiscontinuousPolynomialSpace(grid, order, domains, comp_key)
    elif kind == "P":
        return ContinuousPolynomialSpace(
            grid,
            order,
            domains,
            closed,
            strictly_on_segment,
            element_on_segment,
            comp_key)
    elif kind == "B-P":
        return BarycentricContinuousPolynomialSpace(grid, order, domains, closed, strictly_on_segment, comp_key)
    elif kind == "DUAL":
        if order == 0:
            return DualConstantSpace(grid, domains, closed, strictly_on_segment, comp_key)
        elif order == 1:
            return DualLinearSpace(grid, domains, closed, strictly_on_segment, comp_key)
        else:
            raise ValueError("Only order zero and one dual spaces supported.")
    elif kind == "RT":
        if order > 0:
            raise ValueError(
                "Only order zero Raviart-Thomas spaces supported.")
        return RTSpace(grid, domains, closed, strictly_on_segment, comp_key)
    elif kind == "NC":
        if order > 0:
            raise ValueError(
                "Only order zero Nedelec spaces supported.")
        return NCSpace(grid, domains, closed, strictly_on_segment, comp_key)
    elif kind == "SNC":
        if order > 0:
            raise ValueError(
                "Only order zero Nedelec spaces supported.")
        return SNCSpace(grid, domains, closed, strictly_on_segment, comp_key)
    elif kind == "RWG":
        if order > 0:
            raise ValueError("Only order zero RWG spaces supported.")
        return RWGSpace(grid, domains, closed, strictly_on_segment, comp_key)
    elif kind == "B-RT":
        if order > 0:
            raise ValueError(
                "Only order zero Raviart-Thomas functions supported.")
        return BarycentricRTSpace(grid, domains, closed, strictly_on_segment, comp_key)
    elif kind == "B-NC":
        if order > 0:
            raise ValueError(
                "Only order zero Nedelec functions supported.")
        return BarycentricNCSpace(grid, domains, closed, strictly_on_segment, comp_key)
    elif kind == "B-SNC":
        if order > 0:
            raise ValueError(
                "Only order zero Nedelec functions supported.")
        return BarycentricSNCSpace(grid, domains, closed, strictly_on_segment, comp_key)
    elif kind == "B-RWG":
        if order > 0:
            raise ValueError("Only order zero RWG functions supported.")
        return BarycentricRWGSpace(grid, domains, closed, strictly_on_segment, comp_key)
    elif kind == "BC":
        if order > 0:
            raise ValueError(
                "Only order zero Buffa-Christiansen functions supported.")
        return BuffaChristiansenSpace(grid, domains, closed, strictly_on_segment, comp_key)
    elif kind == "RBC":
        if order > 0:
            raise ValueError(
                "Only order zero Buffa-Christiansen functions supported.")
        return RotatedBuffaChristiansenSpace(grid, domains, closed, strictly_on_segment, comp_key)
    elif kind == "DP-CUSTOM":
        if order == 0 or order == 1:
            if faces_to_include is None:
                faces_to_include = list(range(grid.leaf_view.entity_count(0)))
            if nodes_to_include is None:
                nodes_to_include = list(range(grid.leaf_view.entity_count(2)))
            return CustomDPSpace(grid, order, faces_to_include, nodes_to_include)
        else:
            raise ValueError(
                "Only order zero or one custom DP spaces are supported.")
    else:
        raise ValueError("Unknown space type.")
