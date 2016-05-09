"""Definition of single and double layer Laplace potential operators."""
from . import _common


@_common.potential_logger
def single_layer(space, evaluation_points, parameters=None):
    """Return the Laplace single-layer potential operator

    Parameters
    ----------
    space : bempp.api.space.Space
        The function space over which to assemble the potential.
    evaluation_points : numpy.ndarray
        A (3 x N) array of N evaluation points, where each column corresponds to
        the coordinates of one evaluation point.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.

    """

    import bempp
    from bempp.api.assembly.potential_operator import PotentialOperator
    from bempp.api.assembly.discrete_boundary_operator import GeneralNonlocalDiscreteBoundaryOperator
    from bempp.core.operators.potential.laplace import single_layer_ext

    if space.has_non_barycentric_space:
        space = space.non_barycentric_space

    if parameters is None:
        parameters = bempp.api.global_parameters

    return PotentialOperator(GeneralNonlocalDiscreteBoundaryOperator(single_layer_ext(space._impl, evaluation_points,
                                                                                      parameters)),
                             1, space, evaluation_points)


@_common.potential_logger
def double_layer(space, evaluation_points, parameters=None):
    """Return the Laplace double-layer potential operator

    Parameters
    ----------
    space : bempp.api.space.Space
        The function space over which to assemble the potential.
    evaluation_points : numpy.ndarray
        A (3 x N) array of N evaluation points, where each column corresponds to
        the coordinates of one evaluation point.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.

    """

    import bempp
    from bempp.api.assembly.potential_operator import PotentialOperator
    from bempp.api.assembly.discrete_boundary_operator import GeneralNonlocalDiscreteBoundaryOperator
    from bempp.core.operators.potential.laplace import double_layer_ext

    if space.has_non_barycentric_space:
        space = space.non_barycentric_space

    if parameters is None:
        parameters = bempp.api.global_parameters

    return PotentialOperator(GeneralNonlocalDiscreteBoundaryOperator(double_layer_ext(space._impl, evaluation_points,
                                                                                      parameters)),
                             1, space, evaluation_points)
