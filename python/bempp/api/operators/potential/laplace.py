"""Definition of single and double layer Laplace potential operators."""


def single_layer(space, evaluation_points, parameters=None):
    """Return the Laplace single-layer potential operator."""

    import bempp
    from bempp.api.assembly.potential_operator import PotentialOperator
    from bempp.api.assembly.discrete_boundary_operator import GeneralNonlocalDiscreteBoundaryOperator
    from bempp.core.operators.potential.laplace import single_layer_ext

    if parameters is None:
        parameters = bempp.api.global_parameters

    return PotentialOperator(GeneralNonlocalDiscreteBoundaryOperator(single_layer_ext(space._impl, evaluation_points,
                                                                                      parameters)),
                             1, space, evaluation_points)


def double_layer(space, evaluation_points, parameters=None):
    """Return the Laplace double-layer potential operator."""

    import bempp
    from bempp.api.assembly.potential_operator import PotentialOperator
    from bempp.api.assembly.discrete_boundary_operator import GeneralNonlocalDiscreteBoundaryOperator
    from bempp.core.operators.potential.laplace import double_layer_ext

    if parameters is None:
        parameters = bempp.api.global_parameters

    return PotentialOperator(GeneralNonlocalDiscreteBoundaryOperator(double_layer_ext(space._impl, evaluation_points,
                                                                                      parameters)),
                             1, space, evaluation_points)
