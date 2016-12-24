"""Definition of single and double layer Helmholtz potential operators."""


def single_layer(space, evaluation_points, wave_number, parameters=None):
    """Return the Helmholtz single-layer potential operator

    Parameters
    ----------
    space : bempp.api.space.Space
        The function space over which to assemble the potential.
    evaluation_points : numpy.ndarray
        A (3 x N) array of N evaluation points, where each column corresponds to
        the coordinates of one evaluation point.
    wave_number : complex
        Wavenumber of the operator.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.

    """

    from bempp.api.operators.potential.modified_helmholtz \
        import single_layer as slp

    return slp(
        space, evaluation_points, wave_number / (1j), parameters)


def double_layer(space, evaluation_points, wave_number, parameters=None):
    """Return the Helmholtz double-layer potential operator

    Parameters
    ----------
    space : bempp.api.space.Space
        The function space over which to assemble the potential.
    evaluation_points : numpy.ndarray
        A (3 x N) array of N evaluation points, where each column corresponds to
        the coordinates of one evaluation point.
    wave_number : complex
        Wavenumber of the operator.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given
        the default global parameter object
        `bempp.api.global_parameters` is used.

    """

    from bempp.api.operators.potential.modified_helmholtz \
        import double_layer as dlp

    return dlp(space, evaluation_points, wave_number / (1j), parameters)
