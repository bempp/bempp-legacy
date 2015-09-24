"""Methods for projecting between spaces."""


def discrete_coefficient_projection(space, new_space):
    """Return a discrete operator that projects coefficients from space into new_space."""

    from bempp.api.operators.boundary.sparse import identity
    from bempp.api.assembly.discrete_boundary_operator import InverseSparseDiscreteBoundaryOperator

    ident1 = identity(space, new_space, new_space).weak_form()
    ident2 = identity(new_space, new_space, new_space).weak_form()

    return InverseSparseDiscreteBoundaryOperator(ident2) * ident1

def project_operator(operator, domain=None, range_=None, dual_to_range=None):
    """Project operator onto a new set of spaces."""

    from bempp.api.assembly import _ProjectionBoundaryOperator

    return _ProjectionBoundaryOperator(operator, domain=domain, range_=range, dual_to_range=dual_to_range)

