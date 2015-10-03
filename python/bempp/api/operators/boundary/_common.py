


def multitrace_operator_impl(grid, slp_operator, dlp_operator, hyp_operator, parameters):
    """Implementation of the multitrace operator."""

    import bempp.api
    from bempp.api.assembly import BlockedOperator
    from bempp.api.space import project_operator

    if parameters is None:
        parameters = bempp.api.global_parameters

    blocked_operator = BlockedOperator(2, 2)

    const_space = bempp.api.function_space(grid, "DUAL", 0)
    lin_space = bempp.api.function_space(grid, "P", 1)
    lin_space_bary = bempp.api.function_space(grid, "B-P", 1)
    lin_space_disc_bary = bempp.api.function_space(grid, "B-DP", 1)
    lin_space_disc = bempp.api.function_space(grid.barycentric_grid(), "DP", 1)

    slp = slp_operator(lin_space_disc, lin_space_disc, lin_space_disc, parameters=parameters)
    ident1 = bempp.api.operators.boundary.sparse.identity(lin_space_bary, lin_space_bary, const_space)
    ident2 = bempp.api.operators.boundary.sparse.identity(const_space, const_space, lin_space_bary)
    dlp_disc = dlp_operator(lin_space_disc, lin_space_disc, lin_space_disc, parameters=parameters)
    dlp = project_operator(dlp_disc, domain=lin_space_bary, range_=lin_space_bary,
                           dual_to_range=const_space)

    adlp = project_operator(dlp_disc.transpose(const_space), domain=const_space, range_=const_space,
                            dual_to_range=lin_space_bary)

    blocked_operator[0, 1] = project_operator(slp, domain=const_space, range_=lin_space_bary,
                                              dual_to_range=const_space)

    blocked_operator[1, 0] = (hyp_operator( \
        lin_space_bary, const_space, lin_space_bary,
        use_slp=project_operator(slp, domain=lin_space_disc_bary,
                                 dual_to_range=lin_space_disc_bary),
        parameters=parameters))
    blocked_operator[0, 0] = -dlp
    blocked_operator[1, 1] = adlp

    return blocked_operator
