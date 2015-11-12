# pylint: disable-msg=too-many-arguments

"""Definition of the Maxwell boundary operators."""


def electric_field(domain, range_, dual_to_range,
                   wave_number,
                   label="EFIE", symmetry='no_symmetry',
                   parameters=None, use_slp=False):
    """Return the Maxwell electric field boundary operator.

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
        The electric field operator can be represented as a sparse transformation
        of a Helmholtz single-layer operator. If `use_slp=True` this representation is used.
        It is currently only implemented for the case that the domain, range and dual_to_range
        space are identical. Therefore, the range_ and dual_to_range parameters are ignored.
        If `use_slp=op` for a single-layer boundary operator assembled on a
        suitable space this operator is used to assemble the hypersingular operator.
        Note that if `use_slp=op` is used no checks are performed if the slp operator
        is correctly defined for representing the hypersingular operator. Hence,
        if no care is taken this option can lead to a wrong operator. Also,
        `use_slp=True` or `use_slp=op` is only valid if the `domain` and `dual_to_range`
        spaces are identical.


    """

    import bempp
    from bempp.core.operators.boundary.maxwell import electric_field_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.boundary_operator import BoundaryOperator
    from bempp.api.assembly import LocalBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractLocalOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    if not use_slp:
        return ElementaryBoundaryOperator( \
                ElementaryAbstractIntegralOperator(
            electric_field_ext(parameters, domain._impl, range_._impl, dual_to_range._impl,
                               wave_number, "", symmetry)),
            parameters=parameters, label=label)
    else:

        space = domain

        if not isinstance(use_slp, BoundaryOperator):

            new_space = space.discontinuous_space
            slp = bempp.api.operators.boundary.helmholtz.single_layer(new_space, new_space, new_space, wave_number,
                                                                  parameters=parameters)
        else:
            slp = use_slp

        test_local_ops = []
        trial_local_ops = []

        from bempp.api.assembly.boundary_operator import CompoundBoundaryOperator
        from bempp.core.operators.boundary.sparse import vector_value_times_scalar_ext
        from bempp.core.operators.boundary.sparse import div_times_scalar_ext

        kappa = -1.j * wave_number

        for index in range(3):
            # Definition of range_ does not matter in next operator
            test_local_op = LocalBoundaryOperator(ElementaryAbstractLocalOperator(
                vector_value_times_scalar_ext(slp.dual_to_range._impl, space._impl, space._impl, index)),
                    label='VECTOR_VALUE')
            test_local_ops.append(test_local_op)
            trial_local_ops.append(test_local_op.transpose(space))  # Range parameter arbitrary

        term1 = CompoundBoundaryOperator(test_local_ops, kappa * slp, trial_local_ops, label=label+"_term1")

        test_local_ops = []
        trial_local_ops = []

        div_op = LocalBoundaryOperator(ElementaryAbstractLocalOperator(div_times_scalar_ext(slp.dual_to_range._impl, space._impl, space._impl)),
            label='DIV')
        div_op_transpose = div_op.transpose(space) # Range space does not matter

        term2 = CompoundBoundaryOperator([div_op], (1. / kappa) * slp,
                                         [div_op_transpose], label=label+"_term2")

        return term1 + term2


def magnetic_field(domain, range_, dual_to_range,
                   wave_number,
                   label="MFIE", symmetry='no_symmetry',
                   parameters=None):
    """Return the Maxwell magnetic field boundary operator.

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
    from bempp.core.operators.boundary.maxwell import magnetic_field_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator

    if parameters is None:
        parameters = bempp.api.global_parameters

    return ElementaryBoundaryOperator( \
            ElementaryAbstractIntegralOperator(
        magnetic_field_ext(parameters, domain._impl, range_._impl, dual_to_range._impl,
                           wave_number, "", symmetry)),
        parameters=parameters, label=label)
