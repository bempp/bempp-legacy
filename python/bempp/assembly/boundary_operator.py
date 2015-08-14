"""This module defines various types of boundary operators.

"""


class BoundaryOperator(object):
    """A basic object describing operators acting on boundaries.

    """

    def __init__(self, domain, range_, dual_to_range, parameters=None):

        self._domain = domain
        self._range = range_
        self._dual_to_range = dual_to_range
        self._weak_form = None

        if parameters is None:
            from bempp import global_parameters
            self._parameters = global_parameters
        else:
            self._parameters = parameters

    @property
    def parameters(self):
        """Return the parameters for this boundary operator."""
        return self._parameters

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
        """Return the test space."""
        return self._dual_to_range

    def weak_form(self, recompute=False):
        """Return the discretised weak form."""

        if recompute:
            self._weak_form = None

        if self._weak_form is None:
            self._weak_form = self._weak_form_impl(self._parameters)

        return self._weak_form

    def _weak_form_impl(self, parameters):
        """Returns a weak form. Needs to be implemented by subclasses."""

        raise NotImplementedError("Method not implemented.")

class ElementaryBoundaryOperator(BoundaryOperator):
    """Concrete implementation for elementary integral operators."""

    def __init__(self, abstract_operator, parameters=None):
        super(ElementaryBoundaryOperator, self).__init__(abstract_operator.domain,
                                                         abstract_operator.range,
                                                         abstract_operator.dual_to_range,
                                                         parameters)

        self._impl = abstract_operator

    def _weak_form_impl(self, parameters):
        if parameters.assembly.boundary_operator_assembly_type == 'dense':
            from bempp.assembly.discrete_boundary_operator import \
                DenseDiscreteBoundaryOperator
            return DenseDiscreteBoundaryOperator( \
                    self._impl.assemble_weak_form(parameters).as_matrix())
        else:
            from bempp.assembly.discrete_boundary_operator import \
                    GeneralNonlocalDiscreteBoundaryOperator
            return GeneralNonlocalDiscreteBoundaryOperator( \
                    self._impl.assemble_weak_form(parameters))

class LocalBoundaryOperator(BoundaryOperator):
    """Concrete implementation for local (sparse) boundary operators."""

    def __init__(self, abstract_operator, parameters=None):
        super(LocalBoundaryOperator, self).__init__(abstract_operator.domain,
                                                    abstract_operator.range,
                                                    abstract_operator.dual_to_range,
                                                    parameters)

        self._impl = abstract_operator

    def _weak_form_impl(self, parameters):
        from bempp_ext.assembly.discrete_boundary_operator import convert_to_sparse
        from bempp.assembly.discrete_boundary_operator import SparseDiscreteBoundaryOperator
        return SparseDiscreteBoundaryOperator(\
                convert_to_sparse(self._impl.assemble_weak_form(parameters)))






