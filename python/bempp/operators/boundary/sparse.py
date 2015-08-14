#pylint: disable-msg=too-many-arguments

"""Definition of sparse boundary operators."""


def identity(domain, range_, dual_to_range,
             label='', symmetry='no_symmetry',
             parameters=None):
    """Return the identity operator."""

    import bempp
    from bempp_ext.operators.boundary.sparse import identity_ext
    from bempp.assembly import LocalBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return LocalBoundaryOperator(\
            identity_ext(parameters, domain, range_,
                         dual_to_range, label, symmetry))

def maxwell_identity(domain, range_, dual_to_range,
                     label='', symmetry='no_symmetry',
                     parameters=None):
    """Return the Maxwell identity operator."""

    import bempp
    from bempp_ext.operators.boundary.sparse import maxwell_identity_ext
    from bempp.assembly import LocalBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return LocalBoundaryOperator(\
            maxwell_identity_ext(parameters, domain, range_,
                                 dual_to_range, label, symmetry))

def laplace_beltrami(domain, range_, dual_to_range,
                     label='', symmetry='no_symmetry',
                     parameters=None):
    """Return the Laplace-Beltrami operator."""

    import bempp
    from bempp_ext.operators.boundary.sparse import laplace_beltrami_ext
    from bempp.assembly import LocalBoundaryOperator

    if parameters is None:
        parameters = bempp.global_parameters

    return LocalBoundaryOperator(\
            laplace_beltrami_ext(parameters, domain, range_,
                                 dual_to_range, label, symmetry))


