# pylint: disable-msg=too-many-arguments

"""Definition of the Maxwell boundary operators."""


def electric_field(domain, range_, dual_to_range,
                   wave_number,
                   label='', symmetry='no_symmetry',
                   parameters=None):
    """Return the electric field boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.maxwell import electric_field_ext
    from bempp.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator( \
        electric_field_ext(parameters, domain, range_, dual_to_range,
                         wave_number, label, symmetry),
        parameters=parameters)


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

