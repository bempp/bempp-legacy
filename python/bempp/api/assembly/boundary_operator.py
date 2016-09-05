"""This module defines various types of boundary operators.

"""


def _start_assembly_message(domain, dual_to_range, assembly_type, label):
    """Create a coherent logger message for operator assembly start."""

    return "{2}. START ASSEMBLY. Dim: ({0},{1}). Assembly Type: {3}".format(
        domain.global_dof_count, dual_to_range.global_dof_count, label, assembly_type)


def _end_assembly_message(label, assembly_time):
    """Create a coherent logger message for operator assembly end."""

    return "{1}. FINISHED ASSEMBLY. Time: {0:.2E} sec.".format(assembly_time, label)


def _end_hmat_assembly_message(label, assembly_time, compression_rate, mem_size):
    """Create a coherent logger message for hmat operator assembly end."""

    return "{1}. FINISHED ASSEMBLY. Time: {0:.2E} sec. Mem Size (Mb): {3:.2E}. Compression: {2:.2E}".format(assembly_time, label,
                                                                                                            compression_rate, mem_size)


class BoundaryOperator(object):
    """A basic object describing operators acting on boundaries.

    This class is not supposed to be instantiated directly.

    """

    def __init__(self, domain, range_, dual_to_range, label=""):

        import bempp.api

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
        """Return the discretised weak form.

        Parameters
        ----------
        recompute : bool
            Usually the weak form is cached. If this parameter is set to
            `true` the weak form is recomputed.

        """

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
            if not self.domain.is_compatible(other.space):
                raise ValueError(
                    "Operator domain space does not match GridFunction space.")
            return GridFunction(self.range, projections=self.weak_form() * other.coefficients,
                                dual_space=self.dual_to_range)
        else:
            return NotImplemented

    def __rmul__(self, other):

        import numpy as np

        if np.isscalar(other):
            return _ScaledBoundaryOperator(self, other)
        else:
            return NotImplemented

    def __neg__(self):

        return self.__mul__(-1.0)

    def __sub__(self, other):

        return self.__add__(-other)

    def strong_form(self, recompute=False):
        """Return a discrete operator  that maps into the range space.

        Parameters
        ----------
        recompute : bool
            Usually the strong form is cached. If this parameter is set to
            `true` the strong form is recomputed.
        """
        if recompute is True:
            self._range_map = None

        if self._range_map is None:
           
            # This is the most frequent case and we cache the mass
            # matrix from the space object.
            if self.range == self.dual_to_range:
                self._range_map = self.dual_to_range.inverse_mass_matrix().weak_form()
            else:
                from bempp.api.assembly import InverseSparseDiscreteBoundaryOperator
                from bempp.api.operators.boundary.sparse import identity

                self._range_map = InverseLocalBoundaryOperator(
                    identity(self.range, self.range, self.dual_to_range)).weak_form()

        return self._range_map * self.weak_form(recompute)

    def _weak_form_impl(self):
        """Returns a weak form. Needs to be implemented by subclasses."""

        raise NotImplementedError

    def transpose(self, range_):
        """Return the transpose of a boundary operator.

        Parameters
        ----------
        range_ : bempp.api.space.Space
            The new range space of the transpose. This can not be
            determined automatically.

        """

        return _TransposeBoundaryOperator(self, range_)

    def adjoint(self, range_):
        """Return the adjoint of a boundary operator.

        Parameters
        ----------
        range_ : bempp.api.space.Space
            The new range space of the transpose. This can not be
            determined automatically.

        """
        return _AdjointBoundaryOperator(self, range_)

    def dual_product(self, other):
        """Return the dual product with another boundary operator.

           Notes
           -----
           The dual product B *_D A acts only only on the test functions
           of A  and is defined as 
           self.adjoint(other.range) * other

        """

        from bempp.api.space.projection import rewrite_operator_spaces

        return self.adjoint(other.range) * rewrite_operator_spaces(other,
                other.domain, self.dual_to_range, other.dual_to_range)

class ZeroBoundaryOperator(BoundaryOperator):
    """A boundary operator that represents a zero operator.

    Parameters
    ----------
    domain : bempp.api.space.Space
        Domain space of the operator.
    range_ : bempp.api.space.Space
        Range space of the operator.
    dual_to_range : bempp.api.space.Space
        Dual space to the range space.

    """

    def __init__(self, domain, range_, dual_to_range):
        super(ZeroBoundaryOperator, self).__init__(
            domain, range_, dual_to_range)

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
    """Concrete implementation for elementary integral operators.

    Parameters
    ----------
    abstract_operator : Abstrct boundary operator object
        Various types of abstract operators are defined in
        ``bempp.api.assembly.abstract_boundary_operator``.
    parameters : bempp.api.common.global_parameters
        An optional parameters object (default is None).
    label: string
        An optional operator label (default is "").
    represent_only_singular_part : bool
        When assembled the operator will only contain components for adjacent or 
        overlapping test and trial functions (default false).
        
    Attributes
    ----------
    parameters : bempp.api.common.global_parameters
        Returns the associated parameters object.
    local_assembler : Local assembler object
        Returns the associated local assembler. See also
        ``bempp.api.assembly.assembler``.

    """

    def __init__(self, abstract_operator, parameters=None, label="",
            assemble_only_singular_part=False):
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

        self._assemble_only_singular_part = assemble_only_singular_part

    @property
    def parameters(self):
        """Return the parameters of the operator."""
        return self._parameters

    @property
    def local_assembler(self):
        """Return the local assembler"""

        return self._impl.make_local_assembler(self._parameters)

    def _weak_form_impl(self):
        import time
        import bempp.api

        assembly_mode = self._parameters.assembly.boundary_operator_assembly_type

        bempp.api.LOGGER.info(_start_assembly_message(self.domain,
                                                      self.dual_to_range,
                                                      assembly_mode, self.label))
        start_time = time.time()

        weak_form = self._impl.assemble_weak_form(self._parameters,
                self._assemble_only_singular_part)

        end_time = time.time()

        if assembly_mode == 'hmat' and not self._assemble_only_singular_part:
            from bempp.api.hmat import hmatrix_interface
            mem_size = hmatrix_interface.mem_size(weak_form) / (1.0 * 1024)
            compression_rate = hmatrix_interface.compression_rate(weak_form)
            bempp.api.LOGGER.info(_end_hmat_assembly_message(
                self.label, end_time - start_time, compression_rate, mem_size))
        else:
            bempp.api.LOGGER.info(_end_assembly_message(
                self.label, end_time - start_time))

        return weak_form


class LocalBoundaryOperator(BoundaryOperator):
    """Concrete implementation for local (sparse) boundary operators.

    Parameters
    ----------
    abstract_operator : Abstrct boundary operator object
        Various types of abstract operators are defined in
        ``bempp.api.assembly.abstract_boundary_operator``.
    parameters : bempp.api.common.global_parameters
        An optional parameters object (default is None).
    label: string
        An optional operator label (default is "").

    Attributes
    ----------
    parameters : bempp.api.common.global_parameters
        Returns the associated parameters object.
    local_assembler : Local assembler object
        Returns the associated local assembler. See also
        ``bempp.api.assembly.assember``.

    """

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

    @property
    def local_assembler(self):
        """Return the local assembler"""

        return self._impl.make_local_assembler(self._parameters)

    def _weak_form_impl(self):

        import time
        import bempp.api

        bempp.api.LOGGER.info(_start_assembly_message(self.domain,
                                                      self.dual_to_range,
                                                      'sparse', self.label))
        start_time = time.time()

        weak_form = self._impl.assemble_weak_form(self._parameters)

        end_time = time.time()
        bempp.api.LOGGER.info(_end_assembly_message(
            self.label, end_time - start_time))

        return weak_form


class _ProjectionBoundaryOperator(BoundaryOperator):
    """Define the projection of an operator onto new spaces."""

    def __init__(self, operator, domain, range_, dual_to_range):

        self._operator = operator
        self._domain = domain
        self._range = range_
        self._dual_to_range = dual_to_range

        new_domain = operator.domain if domain is None else domain
        new_range = operator.range if range_ is None else range_
        new_dual_to_range = operator.dual_to_range if dual_to_range is None else dual_to_range

        super(_ProjectionBoundaryOperator, self).__init__(new_domain, new_range, new_dual_to_range,
                                                          'Proj(' + operator.label + ')')

    def _weak_form_impl(self):

        from bempp.api.operators.boundary.sparse import identity
        from .discrete_boundary_operator import InverseSparseDiscreteBoundaryOperator

        projected_weak_form = self._operator.weak_form()

        if self._dual_to_range is not None and self._dual_to_range != self._operator.dual_to_range:
            ident = identity(self._operator.dual_to_range,
                             self._dual_to_range, self._dual_to_range).weak_form()
            projected_weak_form = ident * InverseSparseDiscreteBoundaryOperator(
                identity(self._operator.dual_to_range, self._operator.dual_to_range,
                         self._operator.dual_to_range).weak_form()) * projected_weak_form

        if self._domain is not None and self._domain != self._operator.domain:
            from bempp.api.space.projection import discrete_coefficient_projection

            projected_weak_form *= discrete_coefficient_projection(
                self._domain, self._operator.domain)

        return projected_weak_form


class _ReinterpretSpacesBoundaryOperator(BoundaryOperator):
    """ Return a boundary operator, where the user visible space objects have changed without changing
        the implementation spaces. """

    def __init__(self, operator, domain, range_, dual_to_range):

        super(_ReinterpretSpacesBoundaryOperator, self).__init__(domain, range_, dual_to_range)
        self._operator = operator

    @property
    def local_assembler(self):

        return self._operator.local_assembler

    def _weak_form_impl(self):

        return self._operator.weak_form()


class _SumBoundaryOperator(BoundaryOperator):
    """Return the sum of two boundary operators."""

    def __init__(self, op1, op2):
        if (not op1.domain.is_compatible(op2.domain) or
                not op1.range.is_compatible(op2.range) or
                not op1.dual_to_range.is_compatible(op2.dual_to_range)):
            raise ValueError("Spaces not compatible.")

        super(_SumBoundaryOperator, self).__init__(
            op1.domain, op1.range, op1.dual_to_range, label="(" + op1.label + " + " + op2.label + ")")

        self._op1 = op1
        self._op2 = op2

    def _weak_form_impl(self):
        return self._op1.weak_form() + self._op2.weak_form()


class _ScaledBoundaryOperator(BoundaryOperator):
    """Scale a boundary operator."""

    def __init__(self, op, alpha):
        super(_ScaledBoundaryOperator, self).__init__(
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

        super(_ProductBoundaryOperator, self).__init__(
            op2.domain, op1.range, op1.dual_to_range, label=op1.label + " * " + op2.label)

        self._op1 = op1
        self._op2 = op2

    def _weak_form_impl(self):
        return self._op1.weak_form() * self._op2.strong_form()


class _TransposeBoundaryOperator(BoundaryOperator):
    """Return the transpose of a boundary operator."""

    def __init__(self, op, range_):
        super(_TransposeBoundaryOperator, self).__init__(
            op.dual_to_range, range_, op.domain, label=op.label + ".T")

        self._op = op

    def _weak_form_impl(self):
        return self._op.weak_form().transpose()


class _AdjointBoundaryOperator(BoundaryOperator):
    """Return the adjoint of a boundary operator."""

    def __init__(self, op, range_):
        super(_AdjointBoundaryOperator, self).__init__(
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
            raise ValueError(
                "There must be the same number of test and trial operators.")

        number_of_ops = len(test_local_ops)

        range_ = test_local_ops[0].range
        dual_to_range = test_local_ops[0].dual_to_range
        domain = trial_local_ops[0].domain
        test_domain = test_local_ops[0].domain
        trial_dual_to_range = trial_local_ops[0].dual_to_range

        for i in range(number_of_ops):
            if (not test_local_ops[i].range.is_compatible(range_) or
                    not test_local_ops[i].dual_to_range.is_compatible(dual_to_range) or
                    not trial_local_ops[i].domain.is_compatible(domain) or
                    not test_domain.is_compatible(test_local_ops[i].domain) or
                    not trial_dual_to_range.is_compatible(trial_local_ops[i].dual_to_range)):
                raise ValueError("Incompatible spaces.")

        if parameters is None:
            from bempp.api import global_parameters

            self._parameters = bempp.api.global_parameters
        else:
            self._parameters = parameters

        super(CompoundBoundaryOperator, self).__init__(
            domain, range_, dual_to_range, label=label)

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

        test_inverse = InverseLocalBoundaryOperator(
            identity(self._test_local_ops[0].domain, self._kernel_op.domain,
                     self._kernel_op.dual_to_range,
                     parameters=self._parameters)).weak_form()
        trial_inverse = InverseLocalBoundaryOperator(identity(self._kernel_op.domain, self._kernel_op.domain,
                                                                       self._trial_local_ops[
                                                                           0].dual_to_range,
                                                                       parameters=self._parameters)).weak_form()

        kernel_discrete_op = self._kernel_op.weak_form()
        for i in range(self._number_of_ops):
            discrete_op += ((self._test_local_ops[i].weak_form() * test_inverse) * kernel_discrete_op *
                            (trial_inverse * self._trial_local_ops[i].weak_form()))
        return discrete_op


class RankOneBoundaryOperator(BoundaryOperator):
    """Define a rank one boundary operator.

    This operator is defined as
    op = <w, 1> * <v, 1>, where v are functions in
    `dual_to_range` and w are functions in `domain`.

    Parameters
    ----------
    domain : bempp.api.space.Space
        Domain space.
    range_ : bempp.api.space.Space
        Range space.
    dual_to_range : bempp.api.space.Space
        Dual space to the range space.
    label : string
        Label for the operator.
    symmetry : string
        Symmetry mode. Possible values are: 'no_symmetry',
        'symmetric', 'hermitian'.
    parameters : bempp.api.common.ParameterList
        Parameters for the operator. If none given the
        default global parameter object `bempp.api.global_parameters`
        is used.

    """

    def __init__(self, domain, range_, dual_to_range, label='',
                 parameters=None):

        import bempp.api

        super(RankOneBoundaryOperator, self).__init__(
            domain, range_, dual_to_range, label=label)

        if parameters is None:
            parameters = bempp.api.global_parameters

        self._domain = domain
        self._range = range_
        self._dual_to_range = dual_to_range
        self._parameters = parameters

    def _weak_form_impl(self):

        import bempp.api

        def one_fun(x, n, domain_index, res):
            res[0] = 1

        one_domain = bempp.api.GridFunction(self._domain, fun=one_fun,
                                            parameters=self._parameters)
        one_dual_to_range = bempp.api.GridFunction(self._dual_to_range, fun=one_fun,
                                                   parameters=self._parameters)

        mass_domain = bempp.api.operators.boundary.sparse.identity(self._domain,
                                                                   self._domain, self._domain)
        mass_dual_to_range = bempp.api.operators.boundary.sparse.identity(self._dual_to_range,
                                                                          self._dual_to_range, self._dual_to_range)

        col = one_dual_to_range.projections(self._dual_to_range)
        row = one_domain.projections(self._domain)

        return bempp.api.assembly.DiscreteRankOneOperator(col, row)

class InverseLocalBoundaryOperator(BoundaryOperator):

    def __init__(self, operator):

        import bempp.api

        if not type(operator) == bempp.api.assembly.boundary_operator.LocalBoundaryOperator:
            raise ValueError("Operator must be of type bempp.api.assembly.boundary_operator.LocalBoundaryOperator.")

        super(InverseLocalBoundaryOperator, self).__init__(
            operator.domain, operator.range, operator.dual_to_range,
            label="INV_SPARSE")

        self._domain = operator.domain
        self._range = operator.range
        self._dual_to_range = operator.dual_to_range
        self._op = operator


    def _weak_form_impl(self):

        from bempp.api import InverseSparseDiscreteBoundaryOperator
        from bempp.api.assembly.discrete_boundary_operator import SparseDiscreteBoundaryOperator
        from bempp.api import as_matrix
        from scipy.sparse import block_diag, csc_matrix
        import numpy as np

        if not (self._op.domain.has_local_support and self._op.dual_to_range.has_local_support):
            return InverseSparseDiscreteBoundaryOperator(self._op.weak_form())

        if self._op.domain != self._op.dual_to_range:
            return InverseSparseDiscreteBoundaryOperator(self._op.weak_form())

        # Operator has a block-diagonal form.
        mat = as_matrix(self._op.weak_form()).sorted_indices()

        if mat.getformat != 'csc':
            mat = mat.tocsc()

        if mat.getformat != 'csc':
            mat = mat.tocsc()

        # Figure out, where each block starts and ends

        space = self._domain
        rows = []
        cols = []
        data = []

        indptr = mat.indptr
        indices = mat.indices

        i = 0
        column_indices = []
        diag_data = []

        while i < mat.shape[1]:

            column_indices.append(i)
            block_size = indptr[i+1] - indptr[i]
            block = mat.data[indptr[i]:indptr[i+block_size]]
            diag_data.append(block.reshape((block_size, block_size), order='F'))
            i += block_size

        diag_data = [np.linalg.inv(block) for block in diag_data]

        for i, block in enumerate(diag_data):
            for m in range(block.shape[0]):
                for n in range(block.shape[1]):
                    rows.append(column_indices[i]+m)
                    cols.append(column_indices[i]+n)
                    data.append(block[m,n])


        inv_mat = csc_matrix((data, (rows, cols)), shape=mat.shape,  dtype=np.float64)
        return SparseDiscreteBoundaryOperator(inv_mat)



