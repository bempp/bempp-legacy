"""Definition of single and double layer Helmholtz far field operators."""
from bempp.api.operators.potential import _common


@_common.potential_logger
def single_layer(space, evaluation_points, wave_number, parameters=None):
    """Return the Helmholtz single-layer far field operator.

    Parameters
    ----------
    space : bempp.api.space.Space
        The function space over which to assemble the potential.
    evaluation_points : numpy.ndarray
        A (3 x N) array of N evaluation points, where each column corresponds to
        the coordinates of one evaluation point. For the far field the points
        should lie on the unit sphere.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.

    """

    import bempp
    from bempp.api.assembly.potential_operator import PotentialOperator
    from bempp.api.assembly.discrete_boundary_operator import GeneralNonlocalDiscreteBoundaryOperator
    from bempp.core.operators.far_field.helmholtz import single_layer_ext

    if parameters is None:
        parameters = bempp.api.global_parameters

    return PotentialOperator(GeneralNonlocalDiscreteBoundaryOperator(single_layer_ext(space._impl, evaluation_points,
                                                                                      wave_number,
                                                                                      parameters)),
                             1, space, evaluation_points)


@_common.potential_logger
def double_layer(space, evaluation_points, wave_number, parameters=None):
    """Return the Helmholtz double-layer far field operator.

    Parameters
    ----------
    space : bempp.api.space.Space
        The function space over which to assemble the potential.
    evaluation_points : numpy.ndarray
        A (3 x N) array of N evaluation points, where each column corresponds to
        the coordinates of one evaluation point. For the far field the points
        should lie on the unit sphere.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.

    """

    import bempp
    from bempp.api.assembly.potential_operator import PotentialOperator
    from bempp.api.assembly.discrete_boundary_operator import GeneralNonlocalDiscreteBoundaryOperator
    from bempp.core.operators.far_field.helmholtz import double_layer_ext

    if parameters is None:
        parameters = bempp.api.global_parameters

    return PotentialOperator(GeneralNonlocalDiscreteBoundaryOperator(double_layer_ext(space._impl, evaluation_points,
                                                                                      wave_number,
                                                                                      parameters)),
                             1, space, evaluation_points)
