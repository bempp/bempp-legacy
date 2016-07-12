# pylint: disable-msg=too-many-arguments

"""Definition of the Maxwell boundary operators."""


def electric_field(domain, range_, dual_to_range,
                   wave_number,
                   label="EFIE", symmetry='no_symmetry',
                   parameters=None):
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
    """
#    use_slp : True/False or boundary operator object
#        The electric field operator can be represented as a sparse transformation
#        of a Helmholtz single-layer operator. If `use_slp=True` this representation is used.
#        It is currently only implemented for the case that the domain, range and dual_to_range
#        space are identical. Therefore, the range_ and dual_to_range parameters are ignored.
#        If `use_slp=op` for a single-layer boundary operator assembled on a
#        suitable space this operator is used to assemble the hypersingular operator.
#        Note that if `use_slp=op` is used no checks are performed if the slp operator
#        is correctly defined for representing the hypersingular operator. Hence,
#        if no care is taken this option can lead to a wrong operator. Also,
#        `use_slp=True` or `use_slp=op` is only valid if the `domain` and `dual_to_range`
#        spaces are identical.
#


    import bempp.api
    from bempp.core.operators.boundary.maxwell import electric_field_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.boundary_operator import BoundaryOperator
    from bempp.api.assembly import LocalBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractLocalOperator
    from bempp.api.operators.boundary._common import update_to_non_barycentric_space

    try:
        orig_dual_to_range = dual_to_range
        dual_to_range = dual_to_range._hdiv_space
    except:
        raise ValueError("The dual_to_range space must be a variant of a Nedelec space.")

    new_domain, new_range_, new_dual_to_range = update_to_non_barycentric_space(
        domain, range_, dual_to_range)

    if parameters is None:
        parameters = bempp.api.global_parameters

    efie_op = ElementaryBoundaryOperator(
        ElementaryAbstractIntegralOperator(
            electric_field_ext(parameters, new_domain._impl,
                               new_range_._impl, new_dual_to_range._impl,
                               wave_number, "", symmetry),
            domain, range_, orig_dual_to_range),
        parameters=parameters, label=label)
    return efie_op
#
#        space = domain
#
#        if not isinstance(use_slp, BoundaryOperator):
#
#            disc_space = space.discontinuous_space
#            slp = bempp.api.operators.boundary.helmholtz.single_layer(
#                disc_space, disc_space, disc_space, wave_number,
#                parameters=parameters)
#        else:
#            slp = use_slp
#
#        test_local_ops = []
#        trial_local_ops = []
#
#        from bempp.api.assembly.boundary_operator import CompoundBoundaryOperator
#        from bempp.core.operators.boundary.sparse import vector_value_times_scalar_ext
#        from bempp.core.operators.boundary.sparse import div_times_scalar_ext
#
#        kappa = -1.j * wave_number
#
#        for index in range(3):
#            # Definition of range_ does not matter in next operator
#            test_local_op = LocalBoundaryOperator(ElementaryAbstractLocalOperator(
#                vector_value_times_scalar_ext(
#                    slp.dual_to_range._impl, space._impl, space._impl, index),
#                slp.dual_to_range, space, space),
#                label='VECTOR_VALUE')
#            test_local_ops.append(test_local_op)
#            # Range parameter arbitrary
#            trial_local_ops.append(test_local_op.transpose(space))
#
#        term1 = CompoundBoundaryOperator(
#            test_local_ops, kappa * slp, trial_local_ops, label=label + "_term1")
#
#        test_local_ops = []
#        trial_local_ops = []
#
#        div_op = LocalBoundaryOperator(ElementaryAbstractLocalOperator(div_times_scalar_ext(slp.dual_to_range._impl, space._impl, space._impl),
#                                                                       slp.dual_to_range, space, space),
#                                       label='DIV')
#        div_op_transpose = div_op.transpose(
#            space)  # Range space does not matter
#
#        term2 = CompoundBoundaryOperator([div_op], (1. / kappa) * slp,
#                                         [div_op_transpose], label=label + "_term2")
#
#        efie_op = term1 + term2
#        return efie_op


def calderon_electric_field(grid, wave_number, parameters=None):
    """Return a pair (E^2, E) of the squared EFIE operator E^2 and E itself"""

    import bempp.api

    class EfieSquared(bempp.api.assembly.BoundaryOperator):

        def __init__(self, grid, wave_number, parameters):
            from bempp.api.assembly import InverseSparseDiscreteBoundaryOperator
            from bempp.api.space import project_operator

            bc_space = bempp.api.function_space(grid, "BC", 0)
            tbc_space = bempp.api.function_space(grid, "TBC", 0)
            rwg_space = bempp.api.function_space(grid, "B-RWG", 0)
            snc_space = bempp.api.function_space(grid, "B-SNC", 0)
            rwg_bary_space = bempp.api.function_space(
                grid.barycentric_grid(), "RWG", 0)
            snc_bary_space = bempp.api.function_space(grid.barycentric_grid(), "SNC", 0)
            super(EfieSquared, self).__init__(rwg_space, rwg_space, tbc_space,
                                              label="EFIE_SQUARED")

            self._efie_fine = electric_field(rwg_bary_space, rwg_bary_space, snc_bary_space, wave_number,
                                             parameters=parameters)
            self._efie = project_operator(
                self._efie_fine, domain=rwg_space, range_=rwg_space, dual_to_range=snc_space)
            self._efie2 = project_operator(
                self._efie_fine, domain=bc_space, range_=rwg_space, dual_to_range=tbc_space)
            self._ident = bempp.api.operators.boundary.sparse.identity(
                bc_space, rwg_space, snc_space)
            self._inv_ident = InverseSparseDiscreteBoundaryOperator(
                self._ident.weak_form())

        def _weak_form_impl(self):

            efie_weak = self._efie.weak_form()
            efie2_weak = self._efie2.weak_form()

            return efie2_weak * self._inv_ident * efie_weak

    op = EfieSquared(grid, wave_number, parameters)
    return op, op._efie2


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
    from bempp.api.operators.boundary._common import update_to_non_barycentric_space

    try:
        orig_dual_to_range = dual_to_range
        dual_to_range = dual_to_range._hdiv_space
    except:
        raise ValueError("The dual_to_range space must be a variant of a Nedelec space.")

    new_domain, new_range_, new_dual_to_range = update_to_non_barycentric_space(
        domain, range_, dual_to_range)

    if parameters is None:
        parameters = bempp.api.global_parameters

    mfie_op = ElementaryBoundaryOperator(
        ElementaryAbstractIntegralOperator(
            magnetic_field_ext(parameters, new_domain._impl, new_range_._impl,
                               new_dual_to_range._impl, wave_number, "", symmetry),
            domain, range_, orig_dual_to_range),
        parameters=parameters, label=label)
    return mfie_op
