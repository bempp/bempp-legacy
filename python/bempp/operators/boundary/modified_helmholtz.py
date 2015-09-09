# pylint: disable-msg=too-many-arguments

"""Definition of the modified Helmholtz boundary operators."""


def single_layer(domain, range_, dual_to_range,
                 wave_number,
                 label="SLP", symmetry='no_symmetry',
                 parameters=None):
    """Return the single-layer boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.modified_helmholtz import single_layer_ext
    from bempp.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator( \
        single_layer_ext(parameters, domain, range_,
                         dual_to_range, wave_number, "", symmetry),
        parameters=parameters, label=label)


def double_layer(domain, range_, dual_to_range,
                 wave_number,
                 label="DLP", symmetry='no_symmetry',
                 parameters=None):
    """Return the double-layer boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.modified_helmholtz import double_layer_ext
    from bempp.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator( \
        double_layer_ext(parameters, domain, range_,
                         dual_to_range, wave_number, "", symmetry),
        parameters=parameters, label=label)


def adjoint_double_layer(domain, range_, dual_to_range,
                         wave_number,
                         label="ADJ_DLP", symmetry='no_symmetry',
                         parameters=None):
    """Return the adjoint double-layer boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.modified_helmholtz import adjoint_double_layer_ext
    from bempp.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator( \
        adjoint_double_layer_ext(parameters, domain, range_,
                                 dual_to_range, wave_number, "", symmetry),
        parameters=parameters, label=label)


def hypersingular(domain, range_, dual_to_range, wave_number,
                  label="HYP", symmetry='no_symmetry',
                  parameters=None, use_slp=False):
    """Return the hypersingular boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.modified_helmholtz import hypersingular_ext
    from bempp.assembly.boundary_operator import BoundaryOperator
    from bempp.assembly import ElementaryBoundaryOperator
    from bempp.assembly import LocalBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    if domain != dual_to_range and use_slp:
        print("Compound assembly based on slp operator requires 'domain' and 'dual_to_range' space to be identical." +
              " Switching to standard assembly.")
        use_slp = False

    if not use_slp:
        return ElementaryBoundaryOperator( \
            hypersingular_ext(parameters, domain, range_,
                              dual_to_range, wave_number, "", symmetry),
            parameters=parameters, label=label)
    else:

        if not isinstance(use_slp, BoundaryOperator):
            new_domain = domain.discontinuous_space
            new_dual_to_range = dual_to_range.discontinuous_space
            slp = single_layer(new_domain, range_, new_dual_to_range, wave_number, parameters=parameters)
        else:
            slp = use_slp

        # Test that the spaces are correct.
        if slp.domain != slp.dual_to_range:
            raise ValueError("'domain' and 'dual_to_range' spaces must be identical for the slp operator.")

        if not slp.domain.is_discontinuous:
            raise ValueError("'domain' space of the slp operator must be a discontinuous " +
                             "space of polynomial order larger 0.")

        if not slp.dual_to_range.is_discontinuous:
            raise ValueError("'dual_to_range' space of the slp operator must be a discontinuous " +
                             "space of polynomial order larger 0.")

        # Now generate the compound operator

        test_local_ops = []
        trial_local_ops = []

        from bempp.assembly.boundary_operator import CompoundBoundaryOperator
        from bempp_ext.operators.boundary.sparse import curl_value_ext
        from bempp_ext.operators.boundary.sparse import value_times_normal_ext

        for index in range(3):
            # Definition of range_ does not matter in next operator
            test_local_op = LocalBoundaryOperator(curl_value_ext(slp.dual_to_range, range_, dual_to_range, index))
            test_local_ops.append(test_local_op)
            trial_local_ops.append(test_local_op.transpose(range_))  # Range parameter arbitrary

        term1 = CompoundBoundaryOperator(test_local_ops, slp, trial_local_ops, label=label + "_term1")

        test_local_ops = []
        trial_local_ops = []

        for index in range(3):
            # Definition of range_ does not matter in next operator
            test_local_op = LocalBoundaryOperator(
                value_times_normal_ext(slp.dual_to_range, range_, dual_to_range, index))
            test_local_ops.append(test_local_op)
            trial_local_ops.append(test_local_op.transpose(range_))  # Range parameter arbitrary

        term2 = (wave_number * wave_number) * CompoundBoundaryOperator(test_local_ops, slp, trial_local_ops,
                                                                       label=label + "_term2")

        return term1 + term2
