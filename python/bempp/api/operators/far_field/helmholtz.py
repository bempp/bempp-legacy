"""Definition of single and double layer Helmholtz far field operators."""


def single_layer(space, evaluation_points, wave_number, parameters=None):
    """Return the modified Helmholtz single-layer far field operator."""

    import bempp
    from bempp.api.assembly.potential_operator import PotentialOperator
    from bempp.api.assembly.discrete_boundary_operator import GeneralNonlocalDiscreteBoundaryOperator
    from bempp.core.operators.far_field.helmholtz import single_layer_ext

    if parameters is None:
        parameters = bempp.api.global_parameters

    return PotentialOperator(GeneralNonlocalDiscreteBoundaryOperator(single_layer_ext(space, evaluation_points,
                                                                                      wave_number,
                                                                                      parameters)),
                             1, space, evaluation_points)


def double_layer(space, evaluation_points, wave_number, parameters=None):
    """Return the modified Helmholtz double-layer far field operator."""

    import bempp
    from bempp.api.assembly.potential_operator import PotentialOperator
    from bempp.api.assembly.discrete_boundary_operator import GeneralNonlocalDiscreteBoundaryOperator
    from bempp.core.operators.far_field.helmholtz import double_layer_ext

    if parameters is None:
        parameters = bempp.api.global_parameters

    return PotentialOperator(GeneralNonlocalDiscreteBoundaryOperator(double_layer_ext(space, evaluation_points,
                                                                                      wave_number,
                                                                                      parameters)),
                             1, space, evaluation_points)
