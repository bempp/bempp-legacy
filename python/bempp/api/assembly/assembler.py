"""Various assemblers for integral operators."""


class LocalOperatorLocalAssembler(object):
    """This assembler evaluates local weak forms for local operators."""

    def __init__(self, impl):

        self._impl = impl

    def evaluate_local_weak_forms(self, element_indices):
        """Return local element matrices on the given element indices."""
        return self._impl.evaluate_local_weak_forms(element_indices)


def assemble_dense_block(operator, rows, cols, domain, dual_to_range, parameters=None):
    """Assemble a dense (sub)-block of an elementary integral operator."""

    import bempp
    from .boundary_operator import ElementaryBoundaryOperator
    from .discrete_boundary_operator import DenseDiscreteBoundaryOperator
    from bempp.core.assembly.assembler import assemble_dense_block_ext

    if not isinstance(operator, ElementaryBoundaryOperator):
        raise TypeError(
            "operator must be of type 'ElementaryBoundaryOperator.")

    if parameters is None:
        parameters = operator.parameters

    return DenseDiscreteBoundaryOperator(assemble_dense_block_ext(rows, cols, domain._impl, dual_to_range._impl,
                                                                  operator.local_assembler,
                                                                  parameters).as_matrix())
