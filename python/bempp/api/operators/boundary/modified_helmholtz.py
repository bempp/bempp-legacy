# pylint: disable-msg=too-many-arguments

"""Definition of the modified Helmholtz boundary operators."""


def single_layer(domain, range_, dual_to_range,
                 wave_number,
                 label="SLP", symmetry='no_symmetry',
                 parameters=None):
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

    """

    import bempp
    from bempp.core.operators.boundary.modified_helmholtz import single_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    return ElementaryBoundaryOperator( \
            ElementaryAbstractIntegralOperator(
        single_layer_ext(parameters, domain._impl, range_._impl,
                         dual_to_range._impl, wave_number, "", symmetry)),
        parameters=parameters, label=label)


def double_layer(domain, range_, dual_to_range,
                 wave_number,
                 label="DLP", symmetry='no_symmetry',
                 parameters=None):
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

    """

    import bempp
    from bempp.core.operators.boundary.modified_helmholtz import double_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    return ElementaryBoundaryOperator( \
            ElementaryAbstractIntegralOperator(
        double_layer_ext(parameters, domain._impl, range_._impl,
                         dual_to_range._impl, wave_number, "", symmetry)),
        parameters=parameters, label=label)


def adjoint_double_layer(domain, range_, dual_to_range,
                         wave_number,
                         label="ADJ_DLP", symmetry='no_symmetry',
                         parameters=None):
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

    """

    import bempp
    from bempp.core.operators.boundary.modified_helmholtz import adjoint_double_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    return ElementaryBoundaryOperator( \
            ElementaryAbstractIntegralOperator(
        adjoint_double_layer_ext(parameters, domain._impl, range_._impl,
                                 dual_to_range._impl, wave_number, "", symmetry)),
        parameters=parameters, label=label)


def hypersingular(domain, range_, dual_to_range, wave_number,
                  label="HYP", symmetry='no_symmetry',
                  parameters=None, use_slp=False):
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



    """

    import bempp
    from bempp.core.operators.boundary.modified_helmholtz import hypersingular_ext
    from bempp.api.assembly.boundary_operator import BoundaryOperator
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly import LocalBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractLocalOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    if domain != dual_to_range and use_slp:
        print("Compound assembly based on slp operator requires 'domain' and 'dual_to_range' space to be identical." +
              " Switching to standard assembly.")
        use_slp = False

    if not use_slp:
        return ElementaryBoundaryOperator( \
                ElementaryAbstractIntegralOperator(
            hypersingular_ext(parameters, domain._impl, range_._impl,
                              dual_to_range._impl, wave_number, "", symmetry)),
            parameters=parameters, label=label)
    else:

        if not isinstance(use_slp, BoundaryOperator):
            new_domain = domain.discontinuous_space
            new_dual_to_range = dual_to_range.discontinuous_space
            slp = single_layer(new_domain, range_, new_dual_to_range, wave_number, parameters=parameters)
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
            test_local_op = LocalBoundaryOperator(ElementaryAbstractLocalOperator(curl_value_ext(slp.dual_to_range._impl, range_._impl, dual_to_range._impl, index)),
                label='CURL')
            test_local_ops.append(test_local_op)
            trial_local_ops.append(test_local_op.transpose(range_))  # Range parameter arbitrary

        term1 = CompoundBoundaryOperator(test_local_ops, slp, trial_local_ops, label=label + "_term1")

        test_local_ops = []
        trial_local_ops = []

        for index in range(3):
            # Definition of range_ does not matter in next operator
            test_local_op = LocalBoundaryOperator(ElementaryAbstractLocalOperator(
                value_times_normal_ext(slp.dual_to_range._impl, range_._impl, dual_to_range._impl, index)),
                    label='VALUE_TIMES_NORMAL')
            test_local_ops.append(test_local_op)
            trial_local_ops.append(test_local_op.transpose(range_))  # Range parameter arbitrary

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
        if operator==hypersingular:
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


