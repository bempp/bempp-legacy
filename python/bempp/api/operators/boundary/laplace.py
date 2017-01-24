"""Definition of the Laplace boundary operators."""

# pylint: disable-msg=too-many-arguments

def _single_layer_impl(
        domain, range_, dual_to_range,
        label, symmetry, parameters, assemble_only_singular_part):
    """ Return the actual Laplace single layer operator. """

    from bempp.core.operators.boundary.laplace import single_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import \
        ElementaryAbstractIntegralOperator

    #pylint: disable=protected-access
    return ElementaryBoundaryOperator(
        ElementaryAbstractIntegralOperator(
            single_layer_ext(parameters, domain._impl, range_._impl,
                             dual_to_range._impl, "", symmetry),
            domain, range_, dual_to_range),
        parameters=parameters, label=label,
        assemble_only_singular_part=assemble_only_singular_part)

def _double_layer_impl(
        domain, range_, dual_to_range,
        label, symmetry, parameters, assemble_only_singular_part):
    """ Return the actual Laplace double layer operator. """

    from bempp.core.operators.boundary.laplace import double_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import \
        ElementaryAbstractIntegralOperator

    #pylint: disable=protected-access
    return ElementaryBoundaryOperator(
        ElementaryAbstractIntegralOperator(
            double_layer_ext(parameters, domain._impl, range_._impl,
                             dual_to_range._impl, "", symmetry),
            domain, range_, dual_to_range),
        parameters=parameters, label=label,
        assemble_only_singular_part=assemble_only_singular_part)

def _adjoint_double_layer_impl(
        domain, range_, dual_to_range,
        label, symmetry, parameters, assemble_only_singular_part):
    """ Return the actual Laplace adjoint double layer operator. """

    from bempp.core.operators.boundary.laplace import adjoint_double_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import \
        ElementaryAbstractIntegralOperator

    #pylint: disable=protected-access
    return ElementaryBoundaryOperator(
        ElementaryAbstractIntegralOperator(
            adjoint_double_layer_ext(parameters, domain._impl, range_._impl,
                                     dual_to_range._impl, "", symmetry),
            domain, range_, dual_to_range),
        parameters=parameters, label=label,
        assemble_only_singular_part=assemble_only_singular_part)

def _hypersingular_impl(
        domain, range_, dual_to_range,
        label, symmetry, parameters, assemble_only_singular_part):
    """ Return the actual Laplace hypersingular operator. """

    from bempp.core.operators.boundary.laplace import hypersingular_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import \
        ElementaryAbstractIntegralOperator

    #pylint: disable=protected-access
    return ElementaryBoundaryOperator(
        ElementaryAbstractIntegralOperator(
            hypersingular_ext(parameters, domain._impl, range_._impl,
                              dual_to_range._impl, "", symmetry),
            domain, range_, dual_to_range),
        parameters=parameters, label=label,
        assemble_only_singular_part=assemble_only_singular_part)

def single_layer(domain, range_, dual_to_range,
                 label="SLP", symmetry='no_symmetry',
                 parameters=None, use_projection_spaces=True,
                 assemble_only_singular_part=False):
    """Return the Laplace single-layer boundary operator.

    Parameters
    ----------
    domain : bempp.api.space.Space
        Domain space.
    range_ : bempp.api.space.Space
        Range space.
    dual_to_range : bempp.api.space.Space
        Dual space to the range space.
    label : string
        Label for the operator.
    symmetry : string
        Symmetry mode. Possible values are: 'no_symmetry',
        'symmetric', 'hermitian'.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given the
        default global parameter object `bempp.api.global_parameters`
        is used.
    use_projection_spaces : bool
        Represent operator by projection from higher dimensional space
        if available. This parameter can speed up fast assembly routines,
        such as H-Matrices or FMM (default true).
    assemble_only_singular_part : bool
        When assembled the operator will only contain components for adjacent or
        overlapping test and trial functions (default false).

    """

    from bempp.api.operators.boundary._common import \
        get_operator_with_space_preprocessing

    return get_operator_with_space_preprocessing(
        _single_layer_impl, domain, range_, dual_to_range, label,
        symmetry, parameters, use_projection_spaces,
        assemble_only_singular_part)


def double_layer(domain, range_, dual_to_range,
                 label="DLP", symmetry='no_symmetry',
                 parameters=None,
                 use_projection_spaces=True,
                 assemble_only_singular_part=False):
    """
    Return the Laplace double-layer boundary operator.

    Parameters
    ----------
    domain : bempp.api.space.Space
        Domain space.
    range_ : bempp.api.space.Space
        Range space.
    dual_to_range : bempp.api.space.Space
        Dual space to the range space.
    label : string
        Label for the operator.
    symmetry : string
        Symmetry mode. Possible values are: 'no_symmetry',
        'symmetric', 'hermitian'.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given the
        default global parameter object `bempp.api.global_parameters`
        is used.
    use_projection_spaces : bool
        Represent operator by projection from higher dimensional space
        if available. This parameter can speed up fast assembly routines,
        such as H-Matrices or FMM (default true).
    assemble_only_singular_part : bool
        When assembled the operator will only contain components for adjacent or
        overlapping test and trial functions (default false).

    """

    from bempp.api.operators.boundary._common import \
        get_operator_with_space_preprocessing

    return get_operator_with_space_preprocessing(
        _double_layer_impl, domain, range_, dual_to_range, label,
        symmetry, parameters, use_projection_spaces,
        assemble_only_singular_part)

def adjoint_double_layer(domain, range_, dual_to_range,
                         label="ADJ_DLP", symmetry='no_symmetry',
                         parameters=None,
                         use_projection_spaces=True,
                         assemble_only_singular_part=False):
    """Return the Laplace adjoint double-layer boundary operator.

    Parameters
    ----------
    domain : bempp.api.space.Space
        Domain space.
    range_ : bempp.api.space.Space
        Range space.
    dual_to_range : bempp.api.space.Space
        Dual space to the range space.
    label : string
        Label for the operator.
    symmetry : string
        Symmetry mode. Possible values are: 'no_symmetry',
        'symmetric', 'hermitian'.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given the
        default global parameter object `bempp.api.global_parameters`
        is used.
    use_projection_spaces : bool
        Represent operator by projection from higher dimensional space
        if available. This parameter can speed up fast assembly routines,
        such as H-Matrices or FMM (default true).
    assemble_only_singular_part : bool
        When assembled the operator will only contain components for adjacent or
        overlapping test and trial functions (default false).

    """
    from bempp.api.operators.boundary._common import \
        get_operator_with_space_preprocessing

    return get_operator_with_space_preprocessing(
        _adjoint_double_layer_impl, domain, range_, dual_to_range, label,
        symmetry, parameters, use_projection_spaces,
        assemble_only_singular_part)


#pylint: disable=too-many-locals
def hypersingular(domain, range_, dual_to_range,
                  label="HYP", symmetry='no_symmetry',
                  parameters=None, use_slp=False,
                  use_projection_spaces=True,
                  assemble_only_singular_part=False):
    """Return the Laplace hypersingular boundary operator.

    Parameters
    ----------
    domain : bempp.api.space.Space
        Domain space.
    range_ : bempp.api.space.Space
        Range space.
    dual_to_range : bempp.api.space.Space
        Dual space to the range space.
    label : string
        Label for the operator.
    symmetry : string
        Symmetry mode. Possible values are: 'no_symmetry',
        'symmetric', 'hermitian'.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given the
        default global parameter object `bempp.api.global_parameters`
        is used.
    use_projection_spaces : bool
        Represent operator by projection from higher dimensional space
        if available. This parameter can speed up fast assembly routines,
        such as H-Matrices or FMM (default true).
    use_slp : True/False or boundary operator object
        The hypersingular operator can be represented as a sparse transformation
        of a single-layer operator. If `use_slp=True` this representation
        is used. If `use_slp=op` for a single-layer boundary operator
        assembled on a suitable space this operator is used to assemble the
        hypersingular operator. Note that if `use_slp=op` is used no checks are
        performed if the slp operator is correctly defined for representing
        the hypersingular operator. Hence, if no care is taken this option can
        lead to a wrong operator. Also, `use_slp=True` or `use_slp=op` is only
        valid if the `domain` and `dual_to_range` spaces are identical.
    assemble_only_singular_part : bool
        When assembled the operator will only contain components for adjacent or
        overlapping test and trial functions (default false).
        Note. This option is only used if `use_slp` is not specified.
    """

    import bempp.api
    from bempp.api.assembly.boundary_operator import BoundaryOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    if domain != dual_to_range and use_slp:
        use_slp = False

    if not use_slp:
        from bempp.api.operators.boundary._common import \
            get_operator_with_space_preprocessing
        return get_operator_with_space_preprocessing(
            _hypersingular_impl, domain, range_, dual_to_range, label,
            symmetry, parameters, use_projection_spaces=use_projection_spaces,
            assemble_only_singular_part=assemble_only_singular_part)
    else:

        if domain != dual_to_range:
            raise ValueError(
                "domain and dual_to_range spaces must be identical " +
                "if use_slp=True")

        if not isinstance(use_slp, BoundaryOperator):
            disc_space = domain.discontinuous_space
            slp = single_layer(disc_space, range_,
                               disc_space, parameters=parameters)
        else:
            slp = use_slp


        from bempp.api.assembly.functors import vector_surface_curl_functor
        from bempp.api.assembly.functors import scalar_function_value_functor
        from bempp.api.assembly.functors import \
            single_component_test_trial_integrand_functor
        from bempp.api.space.projection import rewrite_operator_spaces

        compound_op = bempp.api.ZeroBoundaryOperator(
            domain, slp.range, dual_to_range)

        for i in range(3):
            curl_value_op = \
                bempp.api.operators.boundary.sparse.operator_from_functors(
                    domain, slp.domain, slp.domain,
                    scalar_function_value_functor(),
                    vector_surface_curl_functor(),
                    single_component_test_trial_integrand_functor(0, i),
                    label="CURL_OP[{0}]".format(i),
                    parameters=parameters)
            compound_op += curl_value_op.dual_product(slp) * curl_value_op

        # Now generate the compound operator

        return rewrite_operator_spaces(
            compound_op, domain=domain,
            range_=range_, dual_to_range=dual_to_range)


#pylint: disable=invalid-name
def single_layer_and_hypersingular_pair(
        grid, parameters=None, spaces='linear',
        base_slp=None, return_base_slp=False, stabilization_factor=0):
    """Return a pair of hypersingular and single layer operator.

    This function creates a pair of a single-layer and a hypersingular
    operator, where both operators are instantiated using a common
    base single-layer operator. Hence, only one single-layer operator
    needs to be discretized to obtain both operators on the given
    grid.

    Parameters
    ----------
    grid : bempp.api.grid.Grid
        The underlying grid for the multitrace operator
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.
    spaces: string
        Choose 'linear' to assemble the operator
        with continuous linear function spaces for the
        Dirichlet and Neumann component (default). For
        a dual pairing of a linear space for the Dirichlet
        data and piecewise constant space for the Neumann
        data choose 'dual'.
    base_slp : None
        Specify a base single-layer operator to be used. If
        set to None, a base single-layer operator will be
        instantiated by the function.
    return_base_slp : bool
        If True also return the original large space single layer
        operator from which the hypersingular and slp operator
        are derived. Default is False
    stabilization_factor : double
        If not equal to zero add
        this factor times the rank one operator <w 1><v, 1>
        to the hypersingular,
        where w is in the domain space and v in the dual space
        of the hypersingular operator. This regularizes the
        hypersingular operator.

    Returns
    -------
    A pair (slp, hyp) of a single-layer and hypersingular operator.
    If return_base_slp is true a triplet (slp, hyp, base_slp) is
    returned, where base_slp is the single-layer operator, from
    which slp and hyp are obtained via sparse transformations.

    """

    from bempp.api.operators.boundary import _common
    ops = list(_common.slp_and_hyp_impl(
        grid, single_layer, hypersingular, parameters,
        spaces, base_slp, return_base_slp, laplace=True))
    if stabilization_factor != 0:
        from bempp.api.assembly import RankOneBoundaryOperator
        ops[1] += stabilization_factor * RankOneBoundaryOperator(
            ops[1].domain, ops[1].range, ops[1].dual_to_range)
    return ops


def multitrace_operator(grid, parameters=None, spaces='linear', target=None):
    """Return the Laplace multitrace operator.

    Parameters
    ----------
    grid : bempp.api.grid.Grid
        The underlying grid for the multitrace operator
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.
    spaces: string
        Choose 'linear' to assemble the operator
        with continuous linear function spaces for the
        Dirichlet and Neumann component (default). For
        a dual pairing of a linear space for the Dirichlet
        data and piecewise constant space for the Neumann
        data choose 'dual'.
    target: bempp.api.grid.Grid
        Specifies a target grid. If it is different from
        'grid' then the operator maps from 'grid' to
        'target'.

    """

    from bempp.api.operators.boundary import _common
    return _common.multitrace_operator_impl(
        grid, single_layer, double_layer, hypersingular,
        parameters, spaces, laplace=True, target=target)
