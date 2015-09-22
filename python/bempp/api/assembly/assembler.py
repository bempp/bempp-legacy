"""Various assemblers for integral operators."""


def assemble_dense_block(operator, rows, cols, domain, dual_to_range, parameters=None):
    """Assemble a dense (sub)-block of an elementary integral operator."""

    import bempp
    from .boundary_operator import ElementaryBoundaryOperator
    from .discrete_boundary_operator import DenseDiscreteBoundaryOperator
    from bempp.core.assembly.assembler import assemble_dense_block_ext

    if not isinstance(operator, ElementaryBoundaryOperator):
        raise TypeError("operator must be of type 'ElementaryBoundaryOperator.")

    if parameters is None:
        parameters = operator.parameters

    return DenseDiscreteBoundaryOperator(assemble_dense_block_ext(rows, cols, domain, dual_to_range,
                                                                  operator.local_assembler(parameters),
                                                                  parameters).as_matrix())
