def update_to_non_barycentric_space(domain, range_, dual_to_range):
    """Return non-barycentric spaces if possible. Otherwise, return original spaces."""

    if domain.has_non_barycentric_space and dual_to_range.has_non_barycentric_space:
        return (domain.non_barycentric_space, range_, dual_to_range.non_barycentric_space)
    else:
        return (domain, range_, dual_to_range)

def slp_and_hyp_impl(grid, slp_operator, hyp_operator, parameters, spaces='linear', base_slp=None,
        return_base_slp=False, laplace='False'):

    import bempp.api
    from bempp.api.space import project_operator

    if parameters is None:
        parameters = bempp.api.global_parameters

    if spaces=='linear':

        lin_space = bempp.api.function_space(grid, "P", 1)
        lin_space_disc = bempp.api.function_space(grid, "DP", 1)

        if base_slp:
            slp_disc = base_slp
        else:
            slp_disc = slp_operator(lin_space_disc, lin_space_disc, lin_space_disc, parameters=parameters)

        slp = project_operator(slp_disc, domain=lin_space, range_=lin_space, dual_to_range=lin_space)
        hyp = hyp_operator(lin_space, lin_space, lin_space, use_slp=slp_disc, parameters=parameters)

        if return_base_slp:
            return (slp, hyp, slp_disc)
        else:
            return (slp, hyp)


    elif spaces=='dual':

        const_space = bempp.api.function_space(grid, "DUAL", 0)
        lin_space = bempp.api.function_space(grid, "P", 1)
        lin_space_bary = bempp.api.function_space(grid, "B-P", 1)
        lin_space_disc_bary = bempp.api.function_space(grid, "B-DP", 1)
        lin_space_disc = bempp.api.function_space(grid.barycentric_grid(), "DP", 1)

        if base_slp:
            slp = base_slp
        else:
            if laplace:
                const_space_bary = bempp.api.function_space(grid.barycentric_grid(), "DP", 0)
                slp = slp_operator(const_space_bary, const_space_bary, const_space_bary, parameters=parameters)
            else:
                slp = slp_operator(lin_space_disc, lin_space_disc, lin_space_disc, parameters=parameters)

        slp_projected = project_operator(slp, domain=const_space, range_=lin_space_bary,
                                                  dual_to_range=const_space)

        hyp_projected = (hyp_operator( \
            lin_space_bary, const_space, lin_space_bary,
            use_slp=slp,
            parameters=parameters))

        if return_base_slp:
            return (slp_projected, hyp_projected, slp)
        else:
            return (slp_projected, hyp_projected)



def multitrace_operator_impl(grid, slp_operator, dlp_operator, hyp_operator, parameters, 
        spaces='linear', laplace=False):
    if spaces=='linear':
        return multitrace_operator_impl_lin(grid, slp_operator, dlp_operator, hyp_operator,
                parameters)
    elif spaces=='dual':
        return multitrace_operator_impl_dual(grid, slp_operator, dlp_operator, hyp_operator,
                parameters, laplace)
    else:
        raise ValueError("Unknown value for `spaces`.")

def multitrace_operator_impl_lin(grid, slp_operator, dlp_operator, hyp_operator, parameters):
    """Implementation of the multitrace operators for linear spaces."""

    import bempp.api
    from bempp.api.assembly import BlockedOperator
    from bempp.api.space import project_operator

    if parameters is None:
        parameters = bempp.api.global_parameters

    blocked_operator = BlockedOperator(2, 2)

    lin_space = bempp.api.function_space(grid, "P", 1)
    lin_space_disc = bempp.api.function_space(grid, "DP", 1)

    slp_disc = slp_operator(lin_space_disc, lin_space_disc, lin_space_disc, parameters=parameters)
    slp = project_operator(slp_disc, domain=lin_space, range_=lin_space, dual_to_range=lin_space)
    dlp = dlp_operator(lin_space, lin_space, lin_space, parameters=parameters)
    adlp = dlp.transpose(lin_space)
    hyp = hyp_operator(lin_space, lin_space, lin_space, use_slp=slp_disc, parameters=parameters)

    blocked_operator[0, 0] = -dlp
    blocked_operator[0, 1] = slp
    blocked_operator[1, 0] = hyp
    blocked_operator[1, 1] = adlp

    return blocked_operator

def multitrace_operator_impl_dual(grid, slp_operator, dlp_operator, hyp_operator, parameters, laplace):
    """Implementation of the multitrace operator for dual linear/constant spaces."""

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

    if laplace:
        const_space_bary = bempp.api.function_space(grid.barycentric_grid(), "DP", 0)
        slp = slp_operator(const_space_bary, const_space_bary, const_space_bary, parameters=parameters)
    else:
        slp = slp_operator(lin_space_disc, lin_space_disc, lin_space_disc, parameters=parameters)

    dlp_disc = dlp_operator(lin_space_disc, lin_space_disc, lin_space_disc, parameters=parameters)
    dlp = project_operator(dlp_disc, domain=lin_space_bary, range_=lin_space_bary,
                           dual_to_range=const_space)

    adlp = project_operator(dlp_disc.transpose(const_space), domain=const_space, range_=const_space,
                            dual_to_range=lin_space_bary)

    blocked_operator[0, 1] = project_operator(slp, domain=const_space, range_=lin_space_bary,
                                              dual_to_range=const_space)

    blocked_operator[1, 0] = (hyp_operator( \
        lin_space_bary, const_space, lin_space_bary,
        use_slp=slp,
        parameters=parameters))
    blocked_operator[0, 0] = -dlp
    blocked_operator[1, 1] = adlp

    return blocked_operator
