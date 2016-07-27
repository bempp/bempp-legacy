# pylint: disable-msg=too-many-arguments

"""Definition of the Maxwell boundary operators."""

def _electric_field_impl(domain, range_, dual_to_range, wave_number,
        label, symmetry, parameters):
    """ Return the actual electric field operator. """

    from bempp.core.operators.boundary.maxwell import electric_field_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator

    return ElementaryBoundaryOperator(
        ElementaryAbstractIntegralOperator(
            electric_field_ext(parameters, domain._impl, range_._impl,
                             dual_to_range._impl, 
                             wave_number, "", symmetry),
            domain, range_, dual_to_range),
        parameters=parameters, label=label)
    
def _magnetic_field_impl(domain, range_, dual_to_range, wave_number,
        label, symmetry, parameters):
    """ Return the actual magnetic field operator. """

    from bempp.core.operators.boundary.maxwell import magnetic_field_ext
    from bempp.api.assembly import ElementaryBoundaryOperator
    from bempp.api.assembly.abstract_boundary_operator import ElementaryAbstractIntegralOperator

    return ElementaryBoundaryOperator(
        ElementaryAbstractIntegralOperator(
            magnetic_field_ext(parameters, domain._impl, range_._impl,
                             dual_to_range._impl, 
                             wave_number, "", symmetry),
            domain, range_, dual_to_range),
        parameters=parameters, label=label)

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

    from bempp.api.operators.boundary._common import get_wave_operator_with_space_preprocessing
    from bempp.api.space import rewrite_operator_spaces

    try:
        hdiv_dual_to_range = dual_to_range._hdiv_space
    except:
        raise ValueError("The dual space must be a valid Nedelec curl-conforming space.")

    return rewrite_operator_spaces(get_wave_operator_with_space_preprocessing(
            _electric_field_impl, domain, range_, hdiv_dual_to_range, 
            wave_number, label, symmetry, parameters),
            domain, range_, dual_to_range)

def calderon_electric_field(grid, wave_number, parameters=None):
    """Return a pair (E^2, E) of the squared EFIE operator E^2 and E itself"""

    import bempp.api

    class EfieSquared(bempp.api.assembly.BoundaryOperator):

        def __init__(self, grid, wave_number, parameters):
            from bempp.api.assembly import InverseSparseDiscreteBoundaryOperator
            from bempp.api.space import project_operator

            bc_space = bempp.api.function_space(grid, "BC", 0)
            rbc_space = bempp.api.function_space(grid, "RBC", 0)
            rwg_space = bempp.api.function_space(grid, "B-RWG", 0)
            snc_space = bempp.api.function_space(grid, "B-SNC", 0)
            rwg_bary_space = bempp.api.function_space(
                grid.barycentric_grid(), "RWG", 0)
            snc_bary_space = bempp.api.function_space(grid.barycentric_grid(), "SNC", 0)
            super(EfieSquared, self).__init__(rwg_space, rwg_space, rbc_space,
                                              label="EFIE_SQUARED")

            self._efie_fine = electric_field(rwg_bary_space, rwg_bary_space, snc_bary_space, wave_number,
                                             parameters=parameters)
            self._efie = project_operator(
                self._efie_fine, domain=rwg_space, range_=rwg_space, dual_to_range=snc_space)
            self._efie2 = project_operator(
                self._efie_fine, domain=bc_space, range_=rwg_space, dual_to_range=rbc_space)
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

    from bempp.api.operators.boundary._common import get_wave_operator_with_space_preprocessing
    from bempp.api.space import rewrite_operator_spaces

    try:
        hdiv_dual_to_range = dual_to_range._hdiv_space
    except:
        raise ValueError("The dual space must be a valid Nedelec curl-conforming space.")

    return rewrite_operator_spaces(get_wave_operator_with_space_preprocessing(
            _magnetic_field_impl, domain, range_, hdiv_dual_to_range, 
            wave_number, label, symmetry, parameters),
            domain, range_, dual_to_range)
