# pylint: disable-msg=too-many-arguments

"""Definition of sparse boundary operators."""

def operator_from_functors(
        domain, range_, dual_to_range,
        test_functor, trial_functor, integrand_functor,
        label='', symmetry='no_symmetry', parameters=None):
    """Define sparse operator from functors."""
    #pylint: disable=no-name-in-module
    import bempp.api
    from bempp.core.assembly.abstract_boundary_operator import \
        abstract_local_operator_from_functors_ext
    from bempp.api.assembly import LocalBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import \
        ElementaryAbstractLocalOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    #pylint: disable=protected-access
    return LocalBoundaryOperator(
        ElementaryAbstractLocalOperator(
            abstract_local_operator_from_functors_ext(
                domain._impl, range_._impl,
                dual_to_range._impl,
                test_functor._impl, trial_functor._impl,
                integrand_functor._impl, label, symmetry),
            domain, range_, dual_to_range),
        parameters=parameters, label=label)


def identity(domain, range_, dual_to_range,
             label="IDENTITY", symmetry='no_symmetry',
             parameters=None):
    """Return the identity operator.

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

    """
    from bempp.api.assembly.functors import simple_test_trial_integrand_functor
    return operator_from_functors(
        domain, range_, dual_to_range,
        dual_to_range.evaluation_functor,
        domain.evaluation_functor,
        simple_test_trial_integrand_functor(),
        label, symmetry, parameters)

def _maxwell_identity(domain, range_, dual_to_range,
                      label="MAXWELL_IDENTITY", symmetry='no_symmetry',
                      parameters=None):
    """Return the Maxwell identity operator.

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

    """
    from bempp.api.assembly.functors import maxwell_test_trial_integrand_functor
    id_op = operator_from_functors(
        domain, range_, dual_to_range,
        dual_to_range.evaluation_functor,
        domain.evaluation_functor,
        maxwell_test_trial_integrand_functor(),
        label, symmetry, parameters)
    return id_op


def laplace_beltrami(domain, range_, dual_to_range,
                     label="LAPLACE_BELTRAMI", symmetry='no_symmetry',
                     parameters=None):
    """Return the Laplace Beltrami operator.

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

    Notes
    -----
    The spaces for this operator must be spaces of continuous functions.

    """
    from bempp.api.assembly.functors import simple_test_trial_integrand_functor
    from bempp.api.assembly.functors import surface_gradient_functor

    return operator_from_functors(
        domain, range_, dual_to_range,
        surface_gradient_functor(),
        surface_gradient_functor(),
        simple_test_trial_integrand_functor(),
        label, symmetry, parameters)

def multitrace_identity(grid, parameters=None, spaces='linear'):
    """Return the multitrace identity operator.

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
        data choose 'dual'. For the proper dual spaces in
        Maxwell problems choose 'maxwell'.

    """
    from bempp.api.assembly import BlockedOperator
    import bempp.api

    blocked_operator = BlockedOperator(2, 2)

    if spaces == 'dual':
        const_space = bempp.api.function_space(grid, "DUAL", 0)
        lin_space = bempp.api.function_space(grid, "B-P", 1)
        blocked_operator[0, 0] = identity(
            lin_space, lin_space, const_space, parameters=parameters)
        blocked_operator[1, 1] = identity(
            const_space, const_space, lin_space, parameters=parameters)
    elif spaces == 'linear':
        space = bempp.api.function_space(grid, "P", 1)
        blocked_operator[0, 0] = identity(
            space, space, space, parameters=parameters)
        blocked_operator[1, 1] = blocked_operator[0, 0]
    elif spaces == 'maxwell':
        rwg_space = bempp.api.function_space(grid, "B-RWG", 0)
        snc_space = bempp.api.function_space(grid, "B-SNC", 0)
        bc_space = bempp.api.function_space(grid, "BC", 0)
        rbc_space = bempp.api.function_space(grid, "RBC", 0)

        blocked_operator[0, 0] = identity(
            rwg_space, rwg_space, rbc_space, parameters=parameters)
        blocked_operator[1, 1] = identity(
            bc_space, bc_space, snc_space, parameters=parameters)
    else:
        raise ValueError(
            "'spaces' must be one of 'dual', 'linear', or 'maxwell'")

    return blocked_operator
