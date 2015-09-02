# pylint: disable-msg=too-many-arguments

"""Definition of the Laplace boundary operators."""


def single_layer(domain, range_, dual_to_range,
                 label='', symmetry='no_symmetry',
                 parameters=None):
    """Return the single-layer boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.laplace import single_layer_ext
    from bempp.assembly import ElementaryBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator( \
        single_layer_ext(parameters, domain, range_,
                         dual_to_range, label, symmetry),
        parameters=parameters)


def double_layer(domain, range_, dual_to_range,
                 label='', symmetry='no_symmetry',
                 parameters=None):
    """Return the double-layer boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.laplace import double_layer_ext

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator( \
        double_layer_ext(parameters, domain, range_,
                         dual_to_range, label, symmetry),
        parameters=parameters)


def adjoint_double_layer(domain, range_, dual_to_range,
                         label='', symmetry='no_symmetry',
                         parameters=None):
    """Return the adjoint double-layer boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.laplace import adjoint_double_layer_ext

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator( \
        adjoint_double_layer_ext(parameters, domain, range_,
                                 dual_to_range, label, symmetry),
        parameters=parameters)


def hypersingular(domain, range_, dual_to_range,
                  label='', symmetry='no_symmetry',
                  parameters=None):
    """Return the adjoint double-layer boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.laplace import hypersingular_ext

    if parameters is None:
        parameters = bempp.global_parameters

    return ElementaryBoundaryOperator( \
        hypersingular_ext(parameters, domain, range_,
                          dual_to_range, label, symmetry),
        parameters=parameters)
