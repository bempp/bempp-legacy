"""Methods for projecting between spaces."""

def discrete_coefficient_projection(space, new_space):
    """Discrete operator that projects coefficients from space to new_space."""

    from bempp.api.operators.boundary.sparse import identity
    from bempp.api.assembly.discrete_boundary_operator import \
        InverseSparseDiscreteBoundaryOperator

    ident1 = identity(space, new_space, new_space).weak_form()
    ident2 = identity(new_space, new_space, new_space).weak_form()

    return InverseSparseDiscreteBoundaryOperator(ident2) * ident1

def rewrite_operator_spaces(operator, domain, range_, dual_to_range):
    """
    Rewrite the user visible spaces of a boundary operator.

    This method returns a new operator with the given spaces without
    changing the actual internal operator.

    """
    from bempp.api.assembly.boundary_operator import \
        _ReinterpretSpacesBoundaryOperator

    return _ReinterpretSpacesBoundaryOperator(
        operator, domain, range_, dual_to_range)

def project_operator(operator, domain=None, range_=None, dual_to_range=None):
    """Project operator onto a new set of spaces."""

    from bempp.api.assembly.boundary_operator import _ProjectionBoundaryOperator

    if domain is None:
        domain = operator.domain
    if range_ is None:
        range_ = operator.range
    if dual_to_range is None:
        dual_to_range = operator.dual_to_range

    if operator.domain == domain and operator.dual_to_range == dual_to_range:
        if operator.range == range_:
            return operator
        else:
            return rewrite_operator_spaces(
                operator, domain, range_, dual_to_range)

    return _ProjectionBoundaryOperator(
        operator, domain=domain, range_=range_, dual_to_range=dual_to_range)
