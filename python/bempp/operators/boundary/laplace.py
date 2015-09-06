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
                  parameters=None, use_slp=False):
    """Return the hypersingular boundary operator."""

    import bempp
    from bempp_ext.operators.boundary.laplace import hypersingular_ext
    from bempp.assembly.boundary_operator import BoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    if not use_slp:
        return ElementaryBoundaryOperator( \
            hypersingular_ext(parameters, domain, range_,
                              dual_to_range, label, symmetry),
            parameters=parameters)
    else:
        if not isinstance(self, BoundaryOperator):
            slp = single_layer(domain, range_, dual_to_range, parameters=parameters)
        else:
            slp = use_slp
        # Now generate the compound operator

        from bempp.assembly.boundary_operator import ZeroBoundaryOperator
        from bempp_ext.operators.boundary.sparse import curl_value_ext

        op = ZeroBoundaryOperator(domain, range_, dual_to_range)

        for i in range(3):
            curl_op = curl_value_ext(domain, range_, dual_to_range)
            op +=

