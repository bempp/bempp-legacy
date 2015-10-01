"""This module defines various types of boundary operators.

"""


def _start_assembly_message(domain, dual_to_range, assembly_type, label):
    """Create a coherent logger message for operator assembly start."""

    return "Operator: {2}. START ASSEMBLY. Dimensions: ({0},{1}). Assembly Type: {3}".format(
        domain.global_dof_count, dual_to_range.global_dof_count, label, assembly_type)


def _end_assembly_message(label, assembly_time):
    """Create a coherent logger message for operator assembly end."""

    return "Operator: {1}. FINISHED ASSEMBLY. Time: {0} seconds".format(assembly_time, label)


class BoundaryOperator(object):
    """A basic object describing operators acting on boundaries.

    """

    def __init__(self, domain, range_, dual_to_range, label=""):

        self._domain = domain
        self._range = range_
        self._dual_to_range = dual_to_range
        self._weak_form = None
        self._domain_map = None
        self._label = label
        self._range_map = None

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

    @property
    def label(self):
        """Return the label of the operator."""
        return self._label

    def weak_form(self, recompute=False):
        """Return the discretised weak form."""

        if recompute:
            self._weak_form = None

        if self._weak_form is None:
            self._weak_form = self._weak_form_impl()

        return self._weak_form

    def __add__(self, other):
        """Return sum of two boundary operators."""

        return _SumBoundaryOperator(self, other)

    def __mul__(self, other):
        """Return product with a scalar, grid function or other operator."""

        import numpy as np
        from bempp.api import GridFunction

        if np.isscalar(other):
            return _ScaledBoundaryOperator(self, other)
        elif isinstance(other, BoundaryOperator):
            return _ProductBoundaryOperator(self, other)
        elif isinstance(other, GridFunction):
            if self.domain != other.space:
                raise ValueError("Operator domain space does not match GridFunction space.")
            return GridFunction(other.space, coefficients=self.strong_form() * other.coefficients)
        else:
            raise NotImplementedError

    def __rmul__(self, other):

        import numpy as np

        if np.isscalar(other):
            return _ScaledBoundaryOperator(self, other)
        else:
            raise NotImplementedError

    def __neg__(self):

        return self.__mul__(-1.0)

    def __sub__(self, other):

        return self.__add__(-other)

    def strong_form(self):
        """Return a discrete operator  that maps into the range space."""

        if self._range_map is None:
            from bempp.api.operators.boundary.sparse import identity
            from bempp.api.assembly import InverseSparseDiscreteBoundaryOperator

            self._range_map = InverseSparseDiscreteBoundaryOperator( \
                identity(self.range, self.range, self.dual_to_range).weak_form())

        return self._range_map * self.weak_form()

    def _weak_form_impl(self):
        """Returns a weak form. Needs to be implemented by subclasses."""

        raise NotImplementedError

    def transpose(self, range_):
        """Return the transpose of a boundary operator."""

        return _TransposeBoundaryOperator(self, range_)

    def adjoint(self, range_):
        """Return the adjoint of a boundary operator."""

        return _AdjointBoundaryOperator(self, range_)


class ZeroBoundaryOperator(BoundaryOperator):
    """A boundary operator that represents a zero operator."""

    def __init__(self, domain, range_, dual_to_range):
        super(ZeroBoundaryOperator, self).__init__(domain, range_, dual_to_range)

    def _weak_form_impl(self):

        from .discrete_boundary_operator import ZeroDiscreteBoundaryOperator

        return ZeroDiscreteBoundaryOperator(self.dual_to_range.global_dof_count,
                                            self.domain.global_dof_count)

    def __iadd__(self, other):
        if (self.domain != other.domain or
                    self.range != other.range or
                    self.dual_to_range != other.dual_to_range):
            raise ValueError("Spaces not compatible.")

        return other

    def __isub__(self, other):
        if (self.domain != other.domain or
                    self.range != other.range or
                    self.dual_to_range != other.dual_to_range):
            raise ValueError("Spaces not compatible.")

        return -other


class ElementaryBoundaryOperator(BoundaryOperator):
    """Concrete implementation for elementary integral operators."""

    def __init__(self, abstract_operator, parameters=None, label=""):
        super(ElementaryBoundaryOperator, self).__init__(abstract_operator.domain,
                                                         abstract_operator.range,
                                                         abstract_operator.dual_to_range,
                                                         label=label)

        if parameters is None:
            from bempp.api import global_parameters

            self._parameters = global_parameters
        else:
            self._parameters = parameters

        self._impl = abstract_operator

    @property
    def parameters(self):
        """Return the parameters of the operator."""
        return self._parameters

    def local_assembler(self, parameters=None):
        """Return the local assembler"""

        if parameters is None:
            parameters = self.parameters

        return self._impl.make_local_assembler(parameters)

    def _weak_form_impl(self):
        import bempp
        import time

        if self._parameters.assembly.boundary_operator_assembly_type == 'dense':
            from bempp.api.assembly.discrete_boundary_operator import \
                DenseDiscreteBoundaryOperator

            bempp.api.LOGGER.info(_start_assembly_message(self.domain, self.dual_to_range, 'dense', self.label))
            start_time = time.time()

            discrete_operator = DenseDiscreteBoundaryOperator( \
                self._impl.assemble_weak_form(self._parameters).as_matrix())

            end_time = time.time()
            bempp.api.LOGGER.info(_end_assembly_message(self.label, end_time - start_time))

        else:
            from bempp.api.assembly.discrete_boundary_operator import \
                GeneralNonlocalDiscreteBoundaryOperator

            bempp.api.LOGGER.info(_start_assembly_message(self.domain, self.dual_to_range, 'hmat', self.label))
            start_time = time.time()

            discrete_operator = GeneralNonlocalDiscreteBoundaryOperator( \
                self._impl.assemble_weak_form(self._parameters))

            end_time = time.time()

            bempp.api.LOGGER.info(_end_assembly_message(self.label, end_time - start_time))

        return discrete_operator


class LocalBoundaryOperator(BoundaryOperator):
    """Concrete implementation for local (sparse) boundary operators."""

    def __init__(self, abstract_operator, parameters=None, label=""):
        super(LocalBoundaryOperator, self).__init__(abstract_operator.domain,
                                                    abstract_operator.range,
                                                    abstract_operator.dual_to_range,
                                                    label=label)

        if parameters is None:
            from bempp.api import global_parameters

            self._parameters = global_parameters
        else:
            self._parameters = parameters

        self._impl = abstract_operator

    @property
    def parameters(self):
        """Return the parameters of the operator."""

        return self._parameters

    def local_assembler(self, parameters=None):
        """Return the local assembler"""

        if parameters is None:
            parameters = self.parameters

        return self._impl.make_local_assembler(parameters)

    def _weak_form_impl(self):
        from bempp.core.assembly.discrete_boundary_operator import convert_to_sparse
        from bempp.api.assembly.discrete_boundary_operator import SparseDiscreteBoundaryOperator

        import bempp
        import time

        bempp.api.LOGGER.info(_start_assembly_message(self.domain, self.dual_to_range, 'sparse', self.label))
        start_time = time.time()

        discrete_operator = SparseDiscreteBoundaryOperator( \
            convert_to_sparse(self._impl.assemble_weak_form(self._parameters)))

        end_time = time.time()

        bempp.api.LOGGER.info(_end_assembly_message(self.label, end_time - start_time))

        return discrete_operator

class _ProjectionBoundaryOperator(BoundaryOperator):
    """Define the projection of an operator onto new spaces."""

    def __init__(self, operator, domain=None, range_=None, dual_to_range=None):

        self._operator = operator
        self._domain = domain
        self._range = range
        self._dual_to_range = dual_to_range

        new_domain = operator.domain if domain is None else domain
        new_range = operator.range if range_ is None else range_
        new_dual_to_range = operator.dual_to_range if dual_to_range is None else dual_to_range

        super(ProjectionBoundaryOperator, self).__init__(new_domain, new_range, new_dual_to_range,
                                                         'Proj(' + operator.label + ')')

    def _weak_form_impl(self):

        from bempp.api.operators.boundary.sparse import identity
        from .discrete_boundary_operator import InverseSparseDiscreteBoundaryOperator

        projected_weak_form = self._operator.weak_form()

        if self._dual_to_range is not None:

            ident = identity(self._dual_to_range, self._dual_to_range).weak_form()
            result_weak_form = ident * InverseSparseDiscreteBoundaryOperator(
                identity(self._operator._dual_to_range, self._dual_to_range).weak_form()) * projected_weak_form

        if self._domain is not None:

            from bempp.api.space.projection import discrete_coefficient_projection

            projected_weak_form *= discrete_coefficient_projection(self._operator.domain, domain)

        return projected_weak_form

class _SumBoundaryOperator(BoundaryOperator):
    """Return the sum of two boundary operators."""

    def __init__(self, op1, op2):
        if (op1.domain != op2.domain or
                    op1.range != op2.range or
                    op1.dual_to_range != op2.dual_to_range):
            raise ValueError("Spaces not compatible.")

        super(_SumBoundaryOperator, self).__init__( \
            op1.domain, op1.range, op1.dual_to_range, label="(" + op1.label + " + " + op2.label + ")")

        self._op1 = op1
        self._op2 = op2

    def _weak_form_impl(self):
        return self._op1.weak_form() + self._op2.weak_form()


class _ScaledBoundaryOperator(BoundaryOperator):
    """Scale a boundary operator."""

    def __init__(self, op, alpha):
        super(_ScaledBoundaryOperator, self).__init__( \
            op.domain, op.range, op.dual_to_range, label=str(alpha) + " * " + op.label)

        self._op = op
        self._alpha = alpha

    def _weak_form_impl(self):
        return self._op.weak_form() * self._alpha


class _ProductBoundaryOperator(BoundaryOperator):
    """Multiply two boundary operators."""

    def __init__(self, op1, op2):
        if not op2.range.is_compatible(op1.domain):
            raise ValueError("Range space of second operator must be compatible to " +
                             "domain space of first operator.")

        super(_ProductBoundaryOperator, self).__init__( \
            op2.domain, op1.range, op1.dual_to_range, label=op1.label + " * " + op2.label)

        self._op1 = op1
        self._op2 = op2

    def _weak_form_impl(self):
        return self._op1.weak_form() * self._op2.strong_form()


class _TransposeBoundaryOperator(BoundaryOperator):
    """Return the transpose of a boundary operator."""

    def __init__(self, op, range_):
        super(_TransposeBoundaryOperator, self).__init__( \
            op.dual_to_range, range_, op.domain, label=op.label + ".T")

        self._op = op

    def _weak_form_impl(self):
        return self._op.weak_form().transpose()


class _AdjointBoundaryOperator(BoundaryOperator):
    """Return the adjoint of a boundary operator."""

    def __init__(self, op, range_):
        super(_AdjointBoundaryOperator, self).__init__( \
            op.dual_to_range, range_, op.domain, label=op.label + ".H")

        self._op = op

    def _weak_form_impl(self):
        return self._op.weak_form().adjoint()


class CompoundBoundaryOperator(BoundaryOperator):
    """Create a compound boundary operator."""

    def __init__(self, test_local_ops, kernel_op, trial_local_ops,
                 parameters=None, label=""):

        import bempp

        if len(test_local_ops) != len(trial_local_ops):
            raise ValueError("There must be the same number of test and trial operators.")

        number_of_ops = len(test_local_ops)

        range_ = test_local_ops[0].range
        dual_to_range = test_local_ops[0].dual_to_range
        domain = trial_local_ops[0].domain
        test_domain = test_local_ops[0].domain
        trial_dual_to_range = trial_local_ops[0].dual_to_range

        for i in range(number_of_ops):
            if (test_local_ops[i].range != range_ or
                        test_local_ops[i].dual_to_range != dual_to_range or
                        trial_local_ops[i].domain != domain or
                        test_domain != test_local_ops[i].domain or
                        trial_dual_to_range != trial_local_ops[i].dual_to_range):
                raise ValueError("Incompatible spaces.")

        if parameters is None:
            from bempp.api import global_parameters

            self._parameters = bempp.api.global_parameters
        else:
            self._parameters = parameters

        super(CompoundBoundaryOperator, self).__init__(domain, range_, dual_to_range, label=label)

        self._test_local_ops = test_local_ops
        self._trial_local_ops = trial_local_ops
        self._kernel_op = kernel_op
        self._number_of_ops = len(test_local_ops)

    def _weak_form_impl(self):

        from bempp.api.operators.boundary.sparse import identity

        from .discrete_boundary_operator import ZeroDiscreteBoundaryOperator
        from .discrete_boundary_operator import InverseSparseDiscreteBoundaryOperator

        discrete_op = ZeroDiscreteBoundaryOperator(self.dual_to_range.global_dof_count,
                                                   self.domain.global_dof_count)

        test_inverse = InverseSparseDiscreteBoundaryOperator(
            identity(self._test_local_ops[0].domain, self._kernel_op.domain,
                     self._kernel_op.dual_to_range,
                     parameters=self._parameters).weak_form())
        trial_inverse = InverseSparseDiscreteBoundaryOperator(identity(self._kernel_op.domain, self._kernel_op.domain,
                                                                       self._trial_local_ops[0].dual_to_range,
                                                                       parameters=self._parameters).weak_form())

        kernel_discrete_op = self._kernel_op.weak_form()
        for i in range(self._number_of_ops):
            discrete_op += (self._test_local_ops[i].weak_form() * test_inverse * kernel_discrete_op *
                            trial_inverse * self._trial_local_ops[i].weak_form())
        return discrete_op
