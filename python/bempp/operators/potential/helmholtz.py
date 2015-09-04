"""Definition of single and double layer Helmholtz potential operators."""

def single_layer(space, evaluation_points, wave_number, parameters=None):
    """Return the Helmholtz single-layer potential operator"""

    from .modified_helmholtz import single_layer

    return single_layer(space, evaluation_points, wave_number/(1j), parameters)

def double_layer(space, evaluation_points, wave_number, parameters=None):
    """Return the Helmholtz double-layer potential operator"""

    from .modified_helmholtz import double_layer

    return double_layer(space, evaluation_points, wave_number/(1j), parameters)
