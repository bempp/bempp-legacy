# pylint: disable-msg=too-many-arguments

"""Definition of the Laplace boundary operators."""


def single_layer(domain, range_, dual_to_range,
                 label="SLP", symmetry='no_symmetry',
                 parameters=None):
    """Return the single-layer boundary operator."""

    import bempp
    from bempp.core.operators.boundary.laplace import single_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    return ElementaryBoundaryOperator( \
        single_layer_ext(parameters, domain, range_,
                         dual_to_range, "", symmetry),
        parameters=parameters, label=label)


def double_layer(domain, range_, dual_to_range,
                 label="DLP", symmetry='no_symmetry',
                 parameters=None):
    """Return the double-layer boundary operator."""

    import bempp
    from bempp.core.operators.boundary.laplace import double_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator


    if parameters is None:
        parameters = bempp.api.global_parameters

    return ElementaryBoundaryOperator( \
        double_layer_ext(parameters, domain, range_,
                         dual_to_range, "", symmetry),
        parameters=parameters, label=label)


def adjoint_double_layer(domain, range_, dual_to_range,
                         label="ADJ_DLP", symmetry='no_symmetry',
                         parameters=None):
    """Return the adjoint double-layer boundary operator."""

    import bempp
    from bempp.core.operators.boundary.laplace import adjoint_double_layer_ext
    from bempp.api.assembly import ElementaryBoundaryOperator


    if parameters is None:
        parameters = bempp.api.global_parameters

    return ElementaryBoundaryOperator( \
        adjoint_double_layer_ext(parameters, domain, range_,
                                 dual_to_range, "", symmetry),
        parameters=parameters, label=label)


def hypersingular(domain, range_, dual_to_range,
                  label="HYP", symmetry='no_symmetry',
                  parameters=None, use_slp=False):
    """Return the hypersingular boundary operator."""

    import bempp
    from bempp.core.operators.boundary.laplace import hypersingular_ext
    from bempp.api.assembly.boundary_operator import BoundaryOperator
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly import LocalBoundaryOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    if domain != dual_to_range and use_slp:
        print("Compound assembly based on slp operator requires 'domain' and 'dual_to_range' space to be identical." +
              " Switching to standard assembly.")
        use_slp = False

    if not use_slp:
        return ElementaryBoundaryOperator( \
            hypersingular_ext(parameters, domain, range_,
                              dual_to_range, label, symmetry),
            parameters=parameters, label="")
    else:
        if not isinstance(use_slp, BoundaryOperator):
            new_domain = domain.discontinuous_space
            new_dual_to_range = dual_to_range.discontinuous_space
            slp = single_layer(new_domain, range_, new_dual_to_range, parameters=parameters)
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

        from bempp.api.assembly.boundary_operator import CompoundBoundaryOperator
        from bempp.core.operators.boundary.sparse import curl_value_ext

        for index in range(3):
            # Definition of range_ does not matter in next operator
            test_local_op = LocalBoundaryOperator(curl_value_ext(slp.dual_to_range, range_, dual_to_range, index),
                label='CURL')
            test_local_ops.append(test_local_op)
            trial_local_ops.append(test_local_op.transpose(range_)) # Range parameter arbitrary

        return CompoundBoundaryOperator(test_local_ops, slp, trial_local_ops, label=label)


