# pylint: disable-msg=too-many-arguments

"""Definition of the Maxwell boundary operators."""


def electric_field(domain, range_, dual_to_range,
                   wave_number,
                   label='', symmetry='no_symmetry',
                   parameters=None, use_slp=False):
    """Return the electric field boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.maxwell import electric_field_ext
    from bempp.assembly import ElementaryBoundaryOperator
    from bempp.assembly.boundary_operator import BoundaryOperator
    from bempp.assembly import LocalBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    if domain != dual_to_range and use_slp:
        print("Compound assembly based on slp operator requires 'domain' and 'dual_to_range' space to be identical." +
              " Switching to standard assembly.")
        use_slp = False

    if not use_slp:
        return ElementaryBoundaryOperator( \
            electric_field_ext(parameters, domain, range_, dual_to_range,
                               wave_number, label, symmetry),
            parameters=parameters)
    else:
        if not isinstance(use_slp, BoundaryOperator):

            new_domain = domain.discontinuous_space
            new_dual_to_range = dual_to_range.discontinuous_space
            slp = bempp.operators.boundary.helmholtz.single_layer(new_domain, range_, new_dual_to_range, wave_number,
                                                                  parameters=parameters)
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

        test_local_ops = []
        trial_local_ops = []

        from bempp.assembly.boundary_operator import CompoundBoundaryOperator
        from bempp_ext.operators.boundary.sparse import vector_value_times_scalar_ext
        from bempp_ext.operators.boundary.sparse import div_times_scalar_ext

        kappa = -1.j * wave_number

        for index in range(3):
            # Definition of range_ does not matter in next operator
            test_local_op = LocalBoundaryOperator(
                vector_value_times_scalar_ext(slp.dual_to_range, range_, dual_to_range, index))
            test_local_ops.append(test_local_op)
            trial_local_ops.append(test_local_op.transpose(range_))  # Range parameter arbitrary

        term1 = CompoundBoundaryOperator(test_local_ops, kappa * slp, trial_local_ops)

        test_local_ops = []
        trial_local_ops = []

        div_op = LocalBoundaryOperator(div_times_scalar_ext(slp.dual_to_range, range_, dual_to_range))
        div_op_transpose = div_op.transpose(range_) # Range space does not matter

        term2 = CompoundBoundaryOperator([div_op], (1. / kappa) * slp,
                                         [div_op_transpose])

        return term1 + term2


def magnetic_field(domain, range_, dual_to_range,
                   wave_number,
                   label='', symmetry='no_symmetry',
                   parameters=None):
    """Return the magnetic field boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.maxwell import magnetic_field_ext
    from bempp.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator( \
        magnetic_field_ext(parameters, domain, range_, dual_to_range,
                           wave_number, label, symmetry),
        parameters=parameters)
