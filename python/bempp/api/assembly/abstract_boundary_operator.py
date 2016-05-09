"""Interface to abstract boundary operators."""


class ElementaryAbstractIntegralOperator(object):
    """An interfact to abstract elementary non-local boundary operators.

    This object provides methods to discretize non-local boundary
    operators.
    """

    def __init__(self, impl, domain, range_, dual_to_range):
        self._impl = impl
        self._domain = domain
        self._range = range_
        self._dual_to_range = dual_to_range

    def make_local_assembler(self, parameters):
        """Create a local assembler object from the abstract operator."""
        return self._impl.make_local_assembler(parameters)

    def assemble_weak_form(self, parameters):
        """Assemble a boundary integral operator and return the weak form."""
        if parameters.assembly.boundary_operator_assembly_type == 'dense':
            from bempp.api.assembly.discrete_boundary_operator import \
                DenseDiscreteBoundaryOperator

            discrete_operator = DenseDiscreteBoundaryOperator(
                self._impl.assemble_weak_form(parameters).as_matrix())

        else:
            from bempp.api.assembly.discrete_boundary_operator import \
                GeneralNonlocalDiscreteBoundaryOperator

            discrete_operator = GeneralNonlocalDiscreteBoundaryOperator(
                self._impl.assemble_weak_form(parameters))

        return discrete_operator

    @property
    def domain(self):
        """Return the domain space."""
        return self._domain

    @property
    def range(self):
        """Return the range space."""
        return self._range

    @property
    def dual_to_range(self):
        """Return the dual_to_range space."""
        return self._dual_to_range


class ElementaryAbstractLocalOperator(object):
    """An interface to abstract elementary local operators."""

    def __init__(self, impl, domain, range_, dual_to_range):
        self._impl = impl
        self._domain = domain
        self._range = range_
        self._dual_to_range = dual_to_range

    def make_local_assembler(self, parameters):
        """Create a local assembler object from the abstract operator."""
        from .assembler import LocalOperatorLocalAssembler
        return LocalOperatorLocalAssembler(self._impl.make_local_assembler(parameters))

    def assemble_weak_form(self, parameters):
        """Assemble the local operator and return the assembled operator."""
        from bempp.core.assembly.discrete_boundary_operator import convert_to_sparse
        from bempp.api.assembly.discrete_boundary_operator import SparseDiscreteBoundaryOperator

        import bempp.api

        discrete_operator = SparseDiscreteBoundaryOperator( \
            convert_to_sparse(self._impl.assemble_weak_form(parameters)))

        return discrete_operator

    @property
    def domain(self):
        """Return the domain space."""
        return self._domain

    @property
    def range(self):
        """Return the range space."""
        return self._range

    @property
    def dual_to_range(self):
        """Return the dual_to_range space."""
        return self._dual_to_range
