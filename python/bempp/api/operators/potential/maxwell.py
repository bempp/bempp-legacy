"""Definition of electric and magnetic Maxwell potential operators."""


def electric_field(space, evaluation_points, wave_number, parameters=None):
    """Return the Maxwell electric field potential operator."""

    import bempp
    from bempp.api.assembly.potential_operator import PotentialOperator
    from bempp.api.assembly.discrete_boundary_operator import GeneralNonlocalDiscreteBoundaryOperator
    from bempp.core.operators.potential.maxwell import electric_field_ext

    if parameters is None:
        parameters = bempp.api.global_parameters

    return PotentialOperator(GeneralNonlocalDiscreteBoundaryOperator(electric_field_ext(space._impl, evaluation_points,
                                                                                        wave_number,
                                                                                        parameters)),
                             3, space, evaluation_points)


def magnetic_field(space, evaluation_points, wave_number, parameters=None):
    """Return the Maxwell magnetic field potential operator."""

    import bempp
    from bempp.api.assembly.potential_operator import PotentialOperator
    from bempp.api.assembly.discrete_boundary_operator import GeneralNonlocalDiscreteBoundaryOperator
    from bempp.core.operators.potential.maxwell import magnetic_field_ext

    if parameters is None:
        parameters = bempp.api.global_parameters

    return PotentialOperator(GeneralNonlocalDiscreteBoundaryOperator(magnetic_field_ext(space._impl, evaluation_points,
                                                                                        wave_number,
                                                                                        parameters)),
                             3, space, evaluation_points)
