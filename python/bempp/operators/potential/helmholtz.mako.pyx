#cython: embedsignature=True

__all__=['single_layer','double_layer']

from bempp.operators.potential import modified_helmholtz


def single_layer(space, evaluation_points, wave_number,
        parameter_list=None):

    return modified_helmholtz.single_layer(space, evaluation_points,
            wave_number/1j,parameter_list)

def double_layer(space, evaluation_points, wave_number,
        parameter_list=None):

    return modified_helmholtz.double_layer(space, evaluation_points,
            wave_number/1j,parameter_list)
