#pylint: disable-msg=too-many-arguments

"""Definition of the modified Helmholtz boundary operators."""


def single_layer(domain, range_, dual_to_range,
                 wave_number,
                 label='', symmetry='no_symmetry',
                 parameters=None):
    """Return the single-layer boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.modified_helmholtz import single_layer_ext
    from bempp.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator( \
            single_layer_ext(parameters, domain, range_,
                             dual_to_range, wave_number, label, symmetry),
            parameters=parameters)


def double_layer(domain, range_, dual_to_range,
                 wave_number,
                 label='', symmetry='no_symmetry',
                 parameters=None):
    """Return the double-layer boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.modified_helmholtz import double_layer_ext
    from bempp.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator( \
            double_layer_ext(parameters, domain, range_,
                             dual_to_range, wave_number, label, symmetry),
            parameters=parameters)

def adjoint_double_layer(domain, range_, dual_to_range,
                         wave_number,
                         label='', symmetry='no_symmetry',
                         parameters=None):
    """Return the adjoint double-layer boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.modified_helmholtz import adjoint_double_layer_ext
    from bempp.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator( \
            adjoint_double_layer_ext(parameters, domain, range_,
                                     dual_to_range, wave_number, label, symmetry),
            parameters=parameters)

def hypersingular(domain, range_, dual_to_range,
                  wave_number,
                  label='', symmetry='no_symmetry',
                  parameters=None):
    """Return the adjoint double-layer boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.modified_helmholtz import hypersingular_ext
    from bempp.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator( \
            hypersingular_ext(parameters, domain, range_,
                              dual_to_range, wave_number, label, symmetry),
            parameters=parameters)
