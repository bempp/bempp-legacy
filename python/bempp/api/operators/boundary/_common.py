"""Common routines for operators of various types."""

def update_to_non_barycentric_space(domain, range_, dual_to_range):
    """Return non-barycentric spaces if possible (otherwise original spaces)"""

    if (domain.has_non_barycentric_space and
            dual_to_range.has_non_barycentric_space):
        return (domain.non_barycentric_space,
                range_, dual_to_range.non_barycentric_space)
    else:
        return (domain, range_, dual_to_range)

#pylint: disable=invalid-name
def check_for_non_barycentric_spaces(domain, dual_to_range):
    """Check if the spaces can be cast to non-barycentric spaces."""

    if not domain.is_barycentric or not dual_to_range.is_barycentric:
        return False

    return ((domain.has_non_barycentric_space) and
            (dual_to_range.has_non_barycentric_space))


#pylint: disable=too-many-arguments
def get_operator_with_space_preprocessing(
        op_fun, domain, range_, dual_to_range, label,
        symmetry, parameters, use_projection_spaces=True,
        assemble_only_singular_part=False):
    """Assemble operator after preprocessing of spaces."""
    import bempp.api
    from bempp.api.space import rewrite_operator_spaces
    from bempp.api.space import project_operator

    if parameters is None:
        parameters = bempp.api.global_parameters

    if check_for_non_barycentric_spaces(domain, dual_to_range):
        return rewrite_operator_spaces(
            get_operator_with_space_preprocessing(
                op_fun, domain.non_barycentric_space, range_,
                dual_to_range.non_barycentric_space,
                label, symmetry, parameters, use_projection_spaces,
                assemble_only_singular_part),
            domain, range_, dual_to_range)

    if not use_projection_spaces:
        return op_fun(domain, range_, dual_to_range,
                      label, symmetry, parameters, assemble_only_singular_part)
    else:

        operator = op_fun(domain.super_space, range_,
                          dual_to_range.super_space, label,
                          symmetry, parameters, assemble_only_singular_part)

        return project_operator(
            operator, domain=domain, range_=range_,
            dual_to_range=dual_to_range)


def get_wave_operator_with_space_preprocessing(
        op_fun, domain, range_, dual_to_range, wave_number, label,
        symmetry, parameters, use_projection_spaces=True,
        assemble_only_singular_part=False):
    """Preprocess a wave operator."""
    import bempp.api
    from bempp.api.space import rewrite_operator_spaces
    from bempp.api.space import project_operator

    if parameters is None:
        parameters = bempp.api.global_parameters

    if check_for_non_barycentric_spaces(domain, dual_to_range):
        return rewrite_operator_spaces(
            get_wave_operator_with_space_preprocessing(
                op_fun, domain.non_barycentric_space, range_,
                dual_to_range.non_barycentric_space, wave_number,
                label, symmetry, parameters, use_projection_spaces,
                assemble_only_singular_part),
            domain, range_, dual_to_range)

    if not use_projection_spaces:
        return op_fun(
            domain, range_, dual_to_range,
            wave_number, label, symmetry, parameters,
            assemble_only_singular_part)
    else:

        op = op_fun(
            domain.super_space, range_,
            dual_to_range.super_space, wave_number, label,
            symmetry, parameters, assemble_only_singular_part)

        return project_operator(
            op, domain=domain, range_=range_, dual_to_range=dual_to_range)

#pylint: disable=too-many-locals
def slp_and_hyp_impl(
        grid, slp_operator, hyp_operator,
        parameters, spaces='linear', base_slp=None,
        return_base_slp=False, laplace='False'):
    """Assemble slp and hyp operator."""
    import bempp.api
    from bempp.api.space import project_operator

    if parameters is None:
        parameters = bempp.api.global_parameters

    if spaces == 'linear':

        lin_space = bempp.api.function_space(grid, "P", 1)
        lin_space_disc = bempp.api.function_space(grid, "DP", 1)

        if base_slp:
            slp_disc = base_slp
        else:
            slp_disc = slp_operator(
                lin_space_disc, lin_space_disc,
                lin_space_disc, parameters=parameters)

        slp = project_operator(slp_disc, domain=lin_space,
                               range_=lin_space, dual_to_range=lin_space)
        hyp = hyp_operator(lin_space, lin_space, lin_space,
                           use_slp=slp_disc, parameters=parameters)

        if return_base_slp:
            return (slp, hyp, slp_disc)
        else:
            return (slp, hyp)

    elif spaces == 'dual':

        const_space = bempp.api.function_space(grid, "DUAL", 0)
        lin_space = bempp.api.function_space(grid, "P", 1)
        lin_space_bary = bempp.api.function_space(grid, "B-P", 1)
        lin_space_disc = bempp.api.function_space(
            grid.barycentric_grid(), "DP", 1)

        if base_slp:
            slp = base_slp
        else:
            if laplace:
                const_space_bary = bempp.api.function_space(
                    grid.barycentric_grid(), "DP", 0)
                slp = slp_operator(
                    const_space_bary, const_space_bary,
                    const_space_bary, parameters=parameters)
            else:
                slp = slp_operator(lin_space_disc, lin_space_disc,
                                   lin_space_disc, parameters=parameters)

        slp_projected = project_operator(
            slp, domain=const_space, range_=lin_space_bary,
            dual_to_range=const_space)

        hyp_projected = (hyp_operator(
            lin_space_bary, const_space, lin_space_bary,
            use_slp=slp,
            parameters=parameters))

        if return_base_slp:
            return (slp_projected, hyp_projected, slp)
        else:
            return (slp_projected, hyp_projected)


def multitrace_operator_impl(
        grid, slp_operator, dlp_operator, hyp_operator, parameters,
        spaces='linear', laplace=False, target=None):
    """Implementation of multitrace operator."""
    if (target is not None) and (target != grid):
        import bempp.api
        blocked = bempp.api.BlockedOperator(2, 2)
        if spaces == 'linear':

            target_lin = bempp.api.function_space(target, "P", 1)
            source_lin = bempp.api.function_space(grid, "P", 1)

            slp = slp_operator(
                source_lin, target_lin, target_lin, parameters=parameters)
            dlp = dlp_operator(
                source_lin, target_lin, target_lin, parameters=parameters)
            adlp = dlp_operator(
                target_lin, source_lin,
                source_lin, parameters=parameters).transpose(target_lin)
            hyp = hyp_operator(
                source_lin, target_lin, target_lin, parameters=parameters)

            blocked[0, 0] = -dlp
            blocked[0, 1] = slp
            blocked[1, 0] = adlp
            blocked[1, 1] = hyp

        elif spaces == 'dual':

            target_dual = bempp.api.function_space(target, "DUAL", 0)
            target_lin = bempp.api.function_space(target, "B-P", 1)
            source_dual = bempp.api.function_space(grid, "DUAL", 0)
            source_lin = bempp.api.function_space(grid, "B-P", 1)

            slp = slp_operator(
                source_dual, target_lin, target_dual, parameters=parameters)
            dlp = dlp_operator(
                source_lin, target_lin, target_dual, parameters=parameters)
            adlp = dlp_operator(
                target_lin, source_lin, source_dual,
                parameters=parameters).transpose(target_dual)
            hyp = hyp_operator(
                source_lin, target_dual, target_lin, parameters=parameters)

            blocked[0, 0] = -dlp
            blocked[0, 1] = slp
            blocked[1, 0] = adlp
            blocked[1, 1] = hyp

        else:
            raise ValueError("spaces must be one of 'linear' or 'dual'")

        return blocked

    if spaces == 'linear':
        return multitrace_operator_impl_lin(
            grid, slp_operator, dlp_operator, hyp_operator, parameters)
    elif spaces == 'dual':
        return multitrace_operator_impl_dual(
            grid, slp_operator, dlp_operator, hyp_operator,
            parameters, laplace)
    else:
        raise ValueError("Unknown value for `spaces`.")


def multitrace_operator_impl_lin(
        grid, slp_operator, dlp_operator, hyp_operator, parameters):
    """Implementation of the multitrace operators for linear spaces."""
    import bempp.api
    from bempp.api.assembly import BlockedOperator
    from bempp.api.space import project_operator

    if parameters is None:
        parameters = bempp.api.global_parameters

    blocked_operator = BlockedOperator(2, 2)

    lin_space = bempp.api.function_space(grid, "P", 1)
    lin_space_disc = bempp.api.function_space(grid, "DP", 1)

    slp_disc = slp_operator(lin_space_disc, lin_space_disc,
                            lin_space_disc, parameters=parameters)
    slp = project_operator(slp_disc, domain=lin_space,
                           range_=lin_space, dual_to_range=lin_space)
    dlp = dlp_operator(lin_space, lin_space, lin_space, parameters=parameters)
    adlp = dlp.transpose(lin_space)
    hyp = hyp_operator(lin_space, lin_space, lin_space,
                       use_slp=slp_disc, parameters=parameters)

    blocked_operator[0, 0] = -dlp
    blocked_operator[0, 1] = slp
    blocked_operator[1, 0] = hyp
    blocked_operator[1, 1] = adlp

    return blocked_operator


def multitrace_operator_impl_dual(
        grid, slp_operator, dlp_operator, hyp_operator, parameters, laplace):
    """Multitrace operator for dual linear/constant spaces."""
    import bempp.api
    from bempp.api.assembly import BlockedOperator
    from bempp.api.space import project_operator

    if parameters is None:
        parameters = bempp.api.global_parameters

    blocked_operator = BlockedOperator(2, 2)

    const_space = bempp.api.function_space(grid, "DUAL", 0)
    lin_space_bary = bempp.api.function_space(grid, "B-P", 1)
    lin_space_disc = bempp.api.function_space(grid.barycentric_grid(), "DP", 1)

    if laplace:
        const_space_bary = bempp.api.function_space(
            grid.barycentric_grid(), "DP", 0)
        slp = slp_operator(const_space_bary, const_space_bary,
                           const_space_bary, parameters=parameters)
    else:
        slp = slp_operator(lin_space_disc, lin_space_disc,
                           lin_space_disc, parameters=parameters)

    dlp_disc = dlp_operator(lin_space_disc, lin_space_disc,
                            lin_space_disc, parameters=parameters)
    dlp = project_operator(dlp_disc, domain=lin_space_bary,
                           range_=lin_space_bary,
                           dual_to_range=const_space)

    adlp = project_operator(dlp_disc.transpose(const_space),
                            domain=const_space, range_=const_space,
                            dual_to_range=lin_space_bary)

    blocked_operator[0, 1] = project_operator(
        slp, domain=const_space, range_=lin_space_bary,
        dual_to_range=const_space)

    blocked_operator[1, 0] = hyp_operator(
        lin_space_bary, const_space, lin_space_bary, use_slp=slp,
        parameters=parameters)
    blocked_operator[0, 0] = -dlp
    blocked_operator[1, 1] = adlp

    return blocked_operator
