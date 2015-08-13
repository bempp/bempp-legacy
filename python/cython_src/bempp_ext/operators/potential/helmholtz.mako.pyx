#cython: embedsignature=True

__all__=['single_layer','double_layer']

from bempp.operators.potential import modified_helmholtz as _modified_helmholtz


def single_layer(space, evaluation_points, wave_number,
        parameters=None):

    return _modified_helmholtz.single_layer(space, evaluation_points,
            wave_number/1j,parameters)

def double_layer(space, evaluation_points, wave_number,
        parameters=None):

    return _modified_helmholtz.double_layer(space, evaluation_points,
            wave_number/1j,parameters)
