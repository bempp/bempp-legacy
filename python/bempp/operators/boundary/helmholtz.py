#pylint: disable-msg=too-many-arguments

"""Definition of the Helmholtz boundary operators."""


def single_layer(domain, range_, dual_to_range,
                 wave_number,
                 label='', symmetry='no_symmetry',
                 parameters=None):
    """Return the single-layer boundary operator."""

    from .modified_helmholtz import single_layer as sl

    return sl(domain, range_, dual_to_range,
              wave_number/(1j), label, symmetry,
              parameters)


def double_layer(domain, range_, dual_to_range,
                 wave_number,
                 label='', symmetry='no_symmetry',
                 parameters=None):
    """Return the double-layer boundary operator."""

    from .modified_helmholtz import double_layer as dl

    return dl(domain, range_, dual_to_range,
              wave_number/(1j), label, symmetry,
              parameters)

def adjoint_double_layer(domain, range_, dual_to_range,
                         wave_number,
                         label='', symmetry='no_symmetry',
                         parameters=None):
    """Return the adjoint double-layer boundary operator."""

    from .modified_helmholtz import adjoint_double_layer as adl

    return adl(domain, range_, dual_to_range,
               wave_number/(1j), label, symmetry,
               parameters)

def hypersingular(domain, range_, dual_to_range,
                  wave_number,
                  label='', symmetry='no_symmetry',
                  parameters=None):
    """Return the hypersingular boundary operator."""

    from .modified_helmholtz import hypersingular as hyp

    return hyp(domain, range_, dual_to_range,
               wave_number/(1j), label, symmetry,
               parameters)


