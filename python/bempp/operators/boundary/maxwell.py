#pylint: disable-msg=too-many-arguments

"""Definition of the Maxwell boundary operators."""


def single_layer(domain, range_, dual_to_range,
                 wave_number,
                 label='', symmetry='no_symmetry',
                 parameters=None):
    """Return the single-layer boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.maxwell import single_layer_ext
    from bempp.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator(\
            single_layer_ext(parameters, domain, range_, dual_to_range,
                             wave_number, label, symmetry))


def double_layer(domain, range_, dual_to_range,
                 wave_number,
                 label='', symmetry='no_symmetry',
                 parameters=None):
    """Return the double-layer boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.maxwell import double_layer_ext
    from bempp.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator(\
            double_layer_ext(parameters, domain, range_, dual_to_range,
                             wave_number, label, symmetry))

