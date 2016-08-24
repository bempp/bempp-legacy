# pylint: disable-msg=too-many-arguments

"""Definition of the modified Helmholtz boundary operators."""


def _single_layer_impl(domain, range_, dual_to_range, wave_number,
        label, symmetry, parameters, assemble_only_singular_part):
    """ Return the actual modified Helmholtz single layer operator. """

    from bempp.core.operators.boundary.modified_helmholtz import single_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator

    return ElementaryBoundaryOperator(
        ElementaryAbstractIntegralOperator(
            single_layer_ext(parameters, domain._impl, range_._impl,
                             dual_to_range._impl, 
                             wave_number, "", symmetry),
            domain, range_, dual_to_range),
        parameters=parameters, label=label,
        assemble_only_singular_part=assemble_only_singular_part)
    
def _double_layer_impl(domain, range_, dual_to_range, wave_number,
        label, symmetry, parameters, assemble_only_singular_part):
    """ Return the actual modified Helmholtz double layer operator. """

    from bempp.core.operators.boundary.modified_helmholtz import double_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator

    return ElementaryBoundaryOperator(
        ElementaryAbstractIntegralOperator(
            double_layer_ext(parameters, domain._impl, range_._impl,
                             dual_to_range._impl, 
                             wave_number, "", symmetry),
            domain, range_, dual_to_range),
        parameters=parameters, label=label, assemble_only_singular_part=assemble_only_singular_part)

def _adjoint_double_layer_impl(domain, range_, dual_to_range, wave_number,
        label, symmetry, parameters, assemble_only_singular_part):
    """ Return the actual modifed Helmholtz adjoint double layer operator. """

    from bempp.core.operators.boundary.modified_helmholtz import adjoint_double_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator

    return ElementaryBoundaryOperator(
        ElementaryAbstractIntegralOperator(
            adjoint_double_layer_ext(parameters, domain._impl, range_._impl,
                             dual_to_range._impl, wave_number,
                             "", symmetry),
            domain, range_, dual_to_range),
        parameters=parameters, label=label,
        assemble_only_singular_part=assemble_only_singular_part)

def _hypersingular_impl(domain, range_, dual_to_range, wave_number,
        label, symmetry, parameters,
        assemble_only_singular_part):
    """ Return the actual modified Helmholtz hypersingular operator. """

    from bempp.core.operators.boundary.modified_helmholtz import hypersingular_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator

    return ElementaryBoundaryOperator(
        ElementaryAbstractIntegralOperator(
            hypersingular_ext(parameters, domain._impl, range_._impl,
                             dual_to_range._impl, 
                             wave_number, "", symmetry),
            domain, range_, dual_to_range),
        parameters=parameters, label=label,
        assemble_only_singular_part=assemble_only_singular_part)

def single_layer(domain, range_, dual_to_range,
                 wave_number,
                 label="SLP", symmetry='no_symmetry',
                 parameters=None, use_projection_spaces=True,
                 assemble_only_singular_part=False):
    """Return the modified Helmholtz single-layer boundary operator.

    Parameters
    ----------
    domain : bempp.api.space.Space
        Domain space.
    range_ : bempp.api.space.Space
        Range space.
    dual_to_range : bempp.api.space.Space
        Dual space to the range space.
    wave_number : complex
        Wavenumber for the Helmholtz problem.
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

    from bempp.api.operators.boundary._common import get_wave_operator_with_space_preprocessing

    return get_wave_operator_with_space_preprocessing(_single_layer_impl, domain, range_, dual_to_range, 
            wave_number, label, symmetry, parameters, use_projection_spaces, assemble_only_singular_part) 


def double_layer(domain, range_, dual_to_range,
                 wave_number,
                 label="DLP", symmetry='no_symmetry',
                 parameters=None, use_projection_spaces=True,
                 assemble_only_singular_part=False):
    """Return the modified Helmholtz double-layer boundary operator.

    Parameters
    ----------
    domain : bempp.api.space.Space
        Domain space.
    range_ : bempp.api.space.Space
        Range space.
    dual_to_range : bempp.api.space.Space
        Dual space to the range space.
    wave_number : complex
        Wavenumber for the Helmholtz problem.
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

    from bempp.api.operators.boundary._common import get_wave_operator_with_space_preprocessing

    return get_wave_operator_with_space_preprocessing(_double_layer_impl, domain, range_, dual_to_range, 
            wave_number, label, symmetry, parameters, use_projection_spaces, assemble_only_singular_part) 


def adjoint_double_layer(domain, range_, dual_to_range,
                         wave_number,
                         label="ADJ_DLP", symmetry='no_symmetry',
                         parameters=None,
                         use_projection_spaces=True,
                         assemble_only_singular_part=False):
    """Return the modified Helmholtz adjoint double-layer boundary operator.

    Parameters
    ----------
    domain : bempp.api.space.Space
        Domain space.
    range_ : bempp.api.space.Space
        Range space.
    dual_to_range : bempp.api.space.Space
        Dual space to the range space.
    wave_number : complex
        Wavenumber for the Helmholtz problem.
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

    from bempp.api.operators.boundary._common import get_wave_operator_with_space_preprocessing

    return get_wave_operator_with_space_preprocessing(_adjoint_double_layer_impl, domain, range_, dual_to_range, 
            wave_number, label, symmetry, parameters, use_projection_spaces, assemble_only_singular_part) 


def hypersingular(domain, range_, dual_to_range, wave_number,
                  label="HYP", symmetry='no_symmetry',
                  parameters=None, use_slp=False,
                  use_projection_spaces=True,
                  assemble_only_singular_part=False):
    """Return the modified Helmholtz hypersingular boundary operator.

    Parameters
    ----------
    domain : bempp.api.space.Space
        Domain space.
    range_ : bempp.api.space.Space
        Range space.
    dual_to_range : bempp.api.space.Space
        Dual space to the range space.
    wave_number : complex
        Wavenumber for the Helmholtz problem.
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
        of a single-layer operator. If `use_slp=True` this representation is used.
        If `use_slp=op` for a single-layer boundary operator assembled on a
        suitable space this operator is used to assemble the hypersingular operator.
        Note that if `use_slp=op` is used no checks are performed if the slp operator
        is correctly defined for representing the hypersingular operator. Hence,
        if no care is taken this option can lead to a wrong operator. Also,
        `use_slp=True` or `use_slp=op` is only valid if the `domain` and `dual_to_range`
        spaces are identical.
    assemble_only_singular_part : bool
        When assembled the operator will only contain components for adjacent or 
        overlapping test and trial functions (default false).
        Note. This option is only used if `use_slp` is not specified.
    """

    import bempp
    from bempp.core.operators.boundary.modified_helmholtz import hypersingular_ext
    from bempp.api.assembly.boundary_operator import BoundaryOperator
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly import LocalBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractLocalOperator
    from bempp.api.operators.boundary._common import update_to_non_barycentric_space

    if parameters is None:
        parameters = bempp.api.global_parameters

    if domain != dual_to_range and use_slp:
        print("Compound assembly based on slp operator requires 'domain' and 'dual_to_range' space to be identical." +
              " Switching to standard assembly.")
        use_slp = False

    if not use_slp:
        from bempp.api.operators.boundary._common import get_wave_operator_with_space_preprocessing

        return get_wave_operator_with_space_preprocessing(_hypersingular_impl, domain, range_, dual_to_range, 
                wave_number, label, symmetry, parameters, use_projection_spaces, assemble_only_singular_part) 
    else:

        if not isinstance(use_slp, BoundaryOperator):
            disc_domain = domain.discontinuous_space
            disc_dual_to_range = dual_to_range.discontinuous_space
            slp = single_layer(
                disc_domain, range_, disc_dual_to_range, wave_number, parameters=parameters)
        else:
            slp = use_slp

        # Now generate the compound operator

        test_local_ops = []
        trial_local_ops = []

        from bempp.api.assembly.boundary_operator import CompoundBoundaryOperator
        from bempp.core.operators.boundary.sparse import curl_value_ext
        from bempp.core.operators.boundary.sparse import value_times_normal_ext

        for index in range(3):
            # Definition of range_ does not matter in next operator
            test_local_op = LocalBoundaryOperator(
                ElementaryAbstractLocalOperator(
                    curl_value_ext(slp.dual_to_range._impl,
                                   range_._impl, dual_to_range._impl, index),
                    slp.dual_to_range, range_, dual_to_range),
                label='CURL')
            test_local_ops.append(test_local_op)
            trial_local_ops.append(test_local_op.transpose(
                range_))  # Range parameter arbitrary

        term1 = CompoundBoundaryOperator(
            test_local_ops, slp, trial_local_ops, label=label + "_term1")

        test_local_ops = []
        trial_local_ops = []

        for index in range(3):
            # Definition of range_ does not matter in next operator
            test_local_op = LocalBoundaryOperator(ElementaryAbstractLocalOperator(
                value_times_normal_ext(
                    slp.dual_to_range._impl, range_._impl, dual_to_range._impl, index),
                slp.dual_to_range, range_, dual_to_range),
                label='VALUE_TIMES_NORMAL')
            test_local_ops.append(test_local_op)
            trial_local_ops.append(test_local_op.transpose(
                range_))  # Range parameter arbitrary

        term2 = (wave_number * wave_number) * CompoundBoundaryOperator(test_local_ops, slp, trial_local_ops,
                                                                       label=label + "_term2")

        return term1 + term2


def multitrace_operator(grid, wave_number, parameters=None, spaces='linear'):
    """Return the modified Helmholtz multitrace operator.

    Parameters
    ----------
    grid : bempp.api.grid.Grid
        The underlying grid for the multitrace operator
    wave_number : complex
        Wavenumber of the operator.
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

    """

    def op(operator):
        if operator == hypersingular:
            def op_impl(domain, range_, dual_to_range, label="HYP", symmetry="no_symmetry",
                        parameters=None, use_slp=False):
                return hypersingular(domain, range_, dual_to_range, wave_number, label, symmetry, parameters,
                                     use_slp)
            return op_impl
        else:
            import inspect
            defaults = inspect.getargspec(operator).defaults

            def op_impl(domain, range_, dual_to_range, label=defaults[0], symmetry=defaults[1],
                        parameters=None):
                return operator(domain, range_, dual_to_range, wave_number, label, symmetry, parameters)
            return op_impl

    from bempp.api.operators.boundary import _common
    return _common.multitrace_operator_impl(grid, op(single_layer), op(double_layer),
                                            op(hypersingular), parameters, spaces)


def single_layer_and_hypersingular_pair(grid, wave_number, parameters=None, spaces='linear', base_slp=None, return_base_slp=False):
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
    wave_number : complex
        Wavenumber of the operator.
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

    Returns
    -------
    A pair (slp, hyp) of a single-layer and hypersingular operator.
    If return_base_slp is true a triplet (slp, hyp, base_slp) is
    returned, where base_slp is the single-layer operator, from
    which slp and hyp are obtained via sparse transformations.

    """

    from bempp.api.operators.boundary import _common

    def op(operator):
        if operator == hypersingular:
            def op_impl(domain, range_, dual_to_range, label="HYP", symmetry="no_symmetry",
                        parameters=None, use_slp=False):
                return hypersingular(domain, range_, dual_to_range, wave_number, label, symmetry, parameters,
                                     use_slp)
            return op_impl
        else:
            import inspect
            defaults = inspect.getargspec(operator).defaults

            def op_impl(domain, range_, dual_to_range, label=defaults[0], symmetry=defaults[1],
                        parameters=None):
                return operator(domain, range_, dual_to_range, wave_number, label, symmetry, parameters)
            return op_impl

    ops = _common.slp_and_hyp_impl(
        grid, op(single_layer), op(hypersingular), parameters, spaces, base_slp, return_base_slp, laplace=False)
    return ops
