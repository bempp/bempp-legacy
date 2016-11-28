import numpy as _np
from .discrete_boundary_operator import DiscreteBoundaryOperator


def _sum(op1, op2):

    if op1 is None:
        return op2
    elif op2 is None:
        return op1
    else:
        return op1 + op2


def _prod(op1, op2):

    if op1 is None or op2 is None:
        return None
    else:
        return op1 * op2

class BlockedOperatorBase(object):

    def __init__(self, m, n):

        self._m = m
        self._n = n
        self._weak_form = None
        self._range_map = None


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

            _range_ops = _np.empty((self.ndims[0], self.ndims[1]), dtype='O')

            for index in range(self.ndims[0]):

                # This is the most frequent case and we cache the mass
                # matrix from the space object.
                if self.range_spaces[index] == self.dual_to_range_spaces[index]:
                    _range_ops[index, index] = self.dual_to_range_spaces[index].inverse_mass_matrix().weak_form()
                else:
                    from bempp.api.assembly import InverseSparseDiscreteBoundaryOperator
                    from bempp.api.operators.boundary.sparse import identity
                    from bempp.api.assembly.boundary_operator import InverseLocalBoundaryOperator

                    _range_ops[index, index] = InverseLocalBoundaryOperator(
                        identity(self.range_spaces[index], self.range_spaces[index], self.dual_to_range_spaces[index])).weak_form()
            
            self._range_map = BlockedDiscreteOperator(_range_ops)

        return self._range_map * self.weak_form(recompute)

    def __getitem__(self, key):
        raise NotImplementedError()

    @property
    def ndims(self):
        return (self._m, self._n)

    @property
    def range_spaces(self):
        raise NotImplementedError()

    @property
    def dual_to_range_spaces(self):
        raise NotImplementedError()

    @property
    def domain_spaces(self):
        raise NotImplementedError()


    def __add__(self, other):

        if not isinstance(other, BlockedOperatorBase):
            return NotImplementedError

        return SumBlockedOperator(self, other)

        if self.ndims != other.ndims:
            return ValueError("Both blocked operators must have the same dimensions.")

        blocked_operator = BlockedOperator(self.ndims[0], self.ndims[1])

        return blocked_operator

    def __neg__(self):
        return self.__mul__(-1)

    def __sub__(self, other):
        return self.__add__(-other)

    def __mul__(self, other):

        import numpy as np
        import collections

        if np.isscalar(other):

            return ScaledBlockedOperator(self, other)

        elif isinstance(other, BlockedOperatorBase):
            
            return ProductBlockedOperator(self, other)

        elif isinstance(other, collections.Iterable):
            from bempp.api.assembly.grid_function import GridFunction
            list_input = list(other)
            if len(list_input) != self.ndims[1]:
                raise ValueError("Length of input list incompatible with block operator dimension.")
            for item in list_input:
                if not isinstance(item, GridFunction):
                    raise ValueError("All items in the input list must be grid functions.")
            weak_op = self.weak_form()
            x = np.zeros(weak_op.shape[1], dtype=weak_op.dtype)
            col_pos = np.hstack([[0], np.cumsum(weak_op.column_dimensions)])
            row_pos = np.hstack([[0], np.cumsum(weak_op.row_dimensions)])
            for index in range(weak_op.ndims[1]):
                x[col_pos[index]:col_pos[index+1]] = list_input[index].coefficients
            res = weak_op * x
            output_list = []

            for index in range(weak_op.ndims[0]):
                output_list.append(GridFunction(
                    self.range_spaces[index], dual_space=self.dual_to_range_spaces[index],
                    projections=res[row_pos[index]:row_pos[index+1]]))
            return output_list

        else:
            return NotImplementedError

    def __rmul__(self, other):

        import numpy as np
        if np.isscalar(other):
            return self.__mul__(other)
        else:
            return NotImplementedError


class BlockedOperator(BlockedOperatorBase):

    def __init__(self, m, n):

        super(BlockedOperator, self).__init__(m, n)

        self._operators = _np.empty((m, n), dtype=_np.object)
        self._rows = m * [False]
        self._cols = n * [False]

        self._dual_to_range_spaces = m * [None]
        self._range_spaces = m * [None]
        self._domain_spaces = n * [None]

    def __getitem__(self, key):
        import bempp.api

        if self._operators[key] is None:
            return bempp.api.ZeroBoundaryOperator(self.domain_spaces[key[1]],
                    self.range_spaces[key[0]], self.dual_to_range_spaces[key[0]])
        else:
            return self._operators[key]

    def __setitem__(self, key, operator):

        row = key[0]
        col = key[1]

        if self.range_spaces[row] is not None:
            if operator.range != self.range_spaces[row]:
                raise ValueError(
                        "Range space not compatible with self.range_spaces[{0}]".format(row))
        else:
            self._range_spaces[row] = operator.range

        if self.dual_to_range_spaces[row] is not None:
            if operator.dual_to_range != self.dual_to_range_spaces[row]:
                raise ValueError(
                        "Dual to range space not compatible with self.dual_to_range_spaces[{0}]".format(row))
        else:
            self._dual_to_range_spaces[row] = operator.dual_to_range

        if self.domain_spaces[col] is not None:
            if operator.domain != self.domain_spaces[col]:
                raise ValueError(
                        "Domain space {0} not compatible with self.domain_spaces[{1}] = {2}".format(
                            type(operator.domain), col, type(self.domain_spaces[col])))
        else:
            self._domain_spaces[col] = operator.domain

        self._operators[key] = operator
        self._rows[row] = True
        self._cols[col] = True

    def _fill_complete(self):

        return (False not in self._cols) and (False not in self._rows)

    def _weak_form_impl(self):

        if not self._fill_complete():
            raise ValueError(
                "Each row and column must have at least one operator")

        ops = _np.empty((self.ndims[0], self.ndims[1]), dtype='O')

        for i in range(self.ndims[0]):
            for j in range(self.ndims[1]):
                if self._operators[i, j] is not None:
                    ops[i, j] = self._operators[i, j].weak_form()

        return BlockedDiscreteOperator(ops)

    def strong_form(self, recompute=False):

        if not self._fill_complete():
            raise ValueError(
                "Each row and column must have at least one operator")

        return super(BlockedOperator, self).strong_form(recompute)

    @property
    def range_spaces(self):
        return tuple(self._range_spaces)

    @property
    def dual_to_range_spaces(self):
        return tuple(self._dual_to_range_spaces)

    @property
    def domain_spaces(self):
        return tuple(self._domain_spaces)


class SumBlockedOperator(BlockedOperatorBase):
    """Represents the sum of two blocked boundary operators."""

    def __init__(self, op1, op2):

        if op1.ndims != op2.ndims:
            raise ValueError("Incompatible dimensions: {0} != {1}.".format(op1.ndims, op2.ndims))

        for index in range(op1.ndims[0]):
            if op1.range_spaces[index] != op2.range_spaces[index]:
                raise ValueError("Range spaces at index {0} are not identical.".format(index))

        for index in range(op1.ndims[0]):
            if op1.dual_to_range_spaces[index] != op2.dual_to_range_spaces[index]:
                raise ValueError("Dual_to_range spaces at index {0} are not identical.".format(index))

        for index in range(op1.ndims[1]):
            if op1.domain_spaces[index] != op2.domain_spaces[index]:
                raise ValueError("Domain spaces at index {0} are not identical.".format(index))

        self._op1 = op1
        self._op2 = op2

        super(SumBlockedOperator, self).__init__(op1.ndims[0], op1.ndims[1])

    def _weak_form_impl(self):

        return self._op1.weak_form() + self._op2.weak_form()

    def __getitem__(self, key):

        return _sum(self._op1[key], self._op2[key])

    @property
    def range_spaces(self):
        return tuple(self._op1.range_spaces)

    @property
    def dual_to_range_spaces(self):
        return tuple(self._op1.dual_to_range_spaces)

    @property
    def domain_spaces(self):
        return tuple(self._op1.domain_spaces)


class ProductBlockedOperator(BlockedOperatorBase):
    """Represents the Product of two blocked boundary operators."""

    def __init__(self, op1, op2):

        if op1.ndims[1] != op2.ndims[0]:
            raise ValueError("Incompatible dimensions: {0} != {1}.".format(op1.ndims[1], op2.ndims[0]))

        for index in range(op1.ndims[1]):
            if op1.domain_spaces[index] != op2.range_spaces[index]:
                raise ValueError("Range and domain space at index {0} not identical.".format(index))

        self._op1 = op1
        self._op2 = op2

        super(ProductBlockedOperator, self).__init__(op1.ndims[0], op2.ndims[1])

    def _weak_form_impl(self):

        return self._op1.weak_form() * self._op2.strong_form()


    def __getitem__(self, key):

        import bempp.api
        i = key[0]
        j = key[1]

        op = bempp.api.ZeroBoundaryOperator(self.domain_spaces[j],
                self.range_spaces[i], self.dual_to_range_spaces[i])

        for k in range(self._op1.ndims[1]):
            op = _sum(op, _prod(self._op1[i, k], self._op2[k, j]))

        return op

    @property
    def range_spaces(self):
        return tuple(self._op1.range_spaces)

    @property
    def dual_to_range_spaces(self):
        return tuple(self._op1.dual_to_range_spaces)

    @property
    def domain_spaces(self):
        return tuple(self._op2.domain_spaces)

class ScaledBlockedOperator(BlockedOperatorBase):
    """Represents the scalar multiplication of a blocked operator."""

    def __init__(self, op, alpha):

        self._op = op
        self._alpha = alpha

        super(ScaledBlockedOperator, self).__init__(op.ndims[0], op.ndims[1])

    def _weak_form_impl(self):

        return self._alpha * self._op.weak_form()

    def __getitem__(self, key):

        import bempp.api

        if self._op[key] is None:
            return bempp.api.ZeroBoundaryOperator(self.domain_spaces[key[1]], self.range_spaces[key[0]],
                    self.dual_to_range_spaces[key[0]])
        else:
            return self._op[key] * self._alpha

    @property
    def range_spaces(self):
        return tuple(self._op.range_spaces)

    @property
    def dual_to_range_spaces(self):
        return tuple(self._op.dual_to_range_spaces)

    @property
    def domain_spaces(self):
        return tuple(self._op.domain_spaces)



class BlockedDiscreteOperatorBase(DiscreteBoundaryOperator):

    def __init__(self, m, n, dtype, shape):

        self._m = m
        self._n = n
        super(BlockedDiscreteOperatorBase, self).__init__(dtype, shape)

    @property
    def ndims(self):
        return (self._m, self._n)

    @property
    def row_dimensions(self):
        
        raise NotImplementedError()

    @property
    def column_dimensions(self):

        raise NotImplementedError()

    def __getitem__(self, key):

        raise NotImplementedError()

    def _adjoint(self):

        raise NotImplementedError()

    def _transpose(self):

        raise NotImplementedError()

    def __mul__(self, other):

        return self.dot(other)

    def __add__(self, other):

        if isinstance(other, BlockedDiscreteOperatorBase):
            return BlockedDiscreteOperatorSum(self, other)
        else:
            return super(BlockedDiscreteOperatorBase, self).__add__(other)

    def dot(self, other):

        from scipy.sparse.linalg.interface import LinearOperator as _LinearOperator

        if isinstance(other, BlockedDiscreteOperatorBase):
            return BlockedDiscreteOperatorProduct(self, other)
        elif _np.isscalar(other):
            return BlockedScaledDiscreteOperator(self, other)
        elif isinstance(other, _LinearOperator):
            return super(BlockedDiscreteOperatorBase, self).dot(other)
        else:
            x = _np.asarray(other)
            if x.ndim == 1 or (x.ndim == 2 and x.shape[1] == 1):
                return self._matvec(x)
            elif x.ndim == 2:
                return self._matmat(x)
            else:
                raise ValueError("Expect a 1d or 2d array or matrix.")

    def _matmat(self, x):

        return self._matvec(x)

    def __rmul__(self, other):

        if _np.isscalar(other):
            return self * other
        else:
            raise ValueError(
                "Cannot multiply operand of type {0} from the left.".format(type(other)))

    def __call__(self, other):

        return self.dot(other)

    def __matmul__(self, other):

        if np.isscalar(other):
            raise ValueError("Scalar operands not allowed. Use '*' instead.")

        return self.dot(other)

    def __neg__(self):

        return -1 * self

    def __sub__(self, other):

        return self.__add__(-other)

    import scipy
    if scipy.__version__ < '0.16.0':
        def adjoint(self):
            return self._adjoint()
        H = property(adjoint)
        def transpose(self):
            return self._transpose()
        T = property(transpose)

    

class BlockedDiscreteOperator(BlockedDiscreteOperatorBase):

    def __init__(self, ops):

        if not isinstance(ops, _np.ndarray):
            ops = _np.array(ops)

        m = ops.shape[0]
        n = ops.shape[1]

        self._operators = _np.empty((m, n), dtype=_np.object)
        self._rows = _np.zeros(m, dtype=int)
        self._cols = _np.zeros(n, dtype=int)

        for i in range(m):
            for j in range(n):
                if ops[i, j] is None:
                    continue
                if self._rows[i] != 0:
                    if ops[i, j].shape[0] != self._rows[i]:
                        raise ValueError("Block row {0} has incompatible operator sizes.".format(i))
                else:
                    self._rows[i] = ops[i, j].shape[0]

                if self._cols[j] != 0:
                    if ops[i, j].shape[1] != self._cols[j]:
                        raise ValueError("Block column {0} has incompatible operator sizes.".format(j))
                else:
                    self._cols[j] = ops[i, j].shape[1]
                self._operators[i, j] = ops[i, j]

        if not self._fill_complete():
            raise ValueError("Each row and column must contain at least one operator.")

        from bempp.api.assembly.discrete_boundary_operator import ZeroDiscreteBoundaryOperator

        for i in range(m):
            for j in range(n):
                if self._operators[i, j] is None:
                    self._operators[i, j] = ZeroDiscreteBoundaryOperator(self._rows[i], self._cols[j])

        shape = (_np.sum(self._rows), _np.sum(self._cols))

        from bempp.api.utils.data_types import combined_type

        dtype = 'float64'
        for obj in self._operators.ravel():
            if obj is not None:
                dtype = combined_type(dtype, obj.dtype)

        super(BlockedDiscreteOperator, self).__init__(ops.shape[0], ops.shape[1], dtype, shape)

    def __getitem__(self, key):

        return self._operators[key]

    def _fill_complete(self):

        if (0 in self._rows) or (0 in self._cols):
            return False

        return True

    def _matvec(self, x):

        from bempp.api.utils.data_types import combined_type

        if x.ndim == 1:
            x_new = _np.expand_dims(x, 1)
            return self.matvec(x_new).ravel()

        if not self._fill_complete():
            raise ValueError("Not all rows or columns contain operators.")

        row_dim = 0
        res = _np.zeros((self.shape[0], x.shape[1]),
                        dtype=combined_type(self.dtype, x.dtype))

        for i in range(self.ndims[0]):
            col_dim = 0
            local_res = res[row_dim:row_dim + self._rows[i], :]
            for j in range(self.ndims[1]):
                local_x = x[col_dim:col_dim + self._cols[j], :]
                if self._operators[i, j] is not None:
                    op_is_complex = _np.iscomplexobj(
                        self._operators[i, j].dtype.type(1))
                    if _np.iscomplexobj(x) and not op_is_complex:
                        local_res[:] += (self._operators[i, j].dot(_np.real(local_x)) +
                                         1j * self._operators[i, j].dot(_np.imag(local_x)))
                    else:
                        local_res[:] += self._operators[i, j].dot(local_x)
                col_dim += self._cols[j]
            row_dim += self._rows[i]
        return res

    def _get_row_dimensions(self):
        return self._rows

    def _get_column_dimensions(self):
        return self._cols

    def _as_matrix(self):
        if not self._fill_complete():
            raise ValueError("Not all rows or columns contain operators.")
        rows = []
        for i in range(self.ndims[0]):
            row = []
            for j in range(self.ndims[1]):
                if self[i, j] is None:
                    row.append(_np.zeros((self._rows[i], self._cols[j])))
                else:
                    row.append(self[i, j].as_matrix())
            rows.append(_np.hstack(row))
        return _np.vstack(rows)

    row_dimensions = property(_get_row_dimensions)
    column_dimensions = property(_get_column_dimensions)

class BlockedDiscreteOperatorSum(BlockedDiscreteOperatorBase):

    def __init__(self, op1, op2):

        if _np.any(op1.row_dimensions != op2.row_dimensions):
            raise ValueError("Incompatible row dimensions. {0} != {1}".format(
                op1.row_dimensions, op2.row_dimensions))

        if _np.any(op1.column_dimensions != op2.column_dimensions):
            raise ValueError("Incompatible column dimensions. {0} != {1}".format(
                op1.column_dimensions, op2.column_dimensions))

        self._op1 = op1
        self._op2 = op2

        from bempp.api.utils.data_types import combined_type

        super(BlockedDiscreteOperatorSum, self).__init__(
                op1.ndims[0], op1.ndims[1], 
                combined_type(op1.dtype, op2.dtype),
                (op1.shape[0], op1.shape[1]))

    def _matvec(self, x):

        return self._op1 * x + self._op2 * x

    def __getitem__(self, key):

        return self._op1[key] + self._op2[key]

    @property
    def row_dimensions(self):

        return self._op1.row_dimensions

    @property
    def column_dimensions(self):

        return self._op1.column_dimensions


class BlockedDiscreteOperatorProduct(BlockedDiscreteOperatorBase):

    def __init__(self, op1, op2):

        if _np.any(op1.column_dimensions != op2.row_dimensions):
            raise ValueError("Incompatible dimensions. {0} != {1}".format(
                op1.column_dimensions, op2.row_dimensions))

        self._op1 = op1
        self._op2 = op2

        from bempp.api.utils.data_types import combined_type

        super(BlockedDiscreteOperatorProduct, self).__init__(
                op1.ndims[0], op2.ndims[1], 
                combined_type(op1.dtype, op2.dtype),
                (op1.shape[0], op2.shape[1]))

    def _matvec(self, x):

        return self._op1 * (self._op2 * x)

    def __getitem__(self, key):

        from bempp.api.assembly.discrete_boundary_operator import ZeroDiscreteBoundaryOperator

        i = key[0]
        j = key[1]

        op = ZeroDiscreteBoundaryOperator(self.row_dimensions[i], self.column_dimensions[j])

        for k in range(self._op1.ndims[1]):
            op += self._op1[i, k] * self._op2[k, j]

        return op

    @property
    def row_dimensions(self):

        return self._op1.row_dimensions

    @property 
    def column_dimensions(self):

        return self._op2.column_dimensions


class BlockedScaledDiscreteOperator(BlockedDiscreteOperatorBase):

    def __init__(self, op, alpha):

        self._op = op
        self._alpha = alpha

        from bempp.api.utils.data_types import combined_type

        if _np.iscomplex(alpha):
            d = _np.dtype('complex128')
        else:
            d = op.dtype

        super(BlockedScaledDiscreteOperator, self).__init__(
                op.ndims[0], op.ndims[1], d,
                (op.shape[0], op.shape[1]))

    def _matvec(self, x):

        return self._alpha * (self._op * x)

    def __getitem__(self, key):

        return self._op[key] * self._alpha

    @property
    def row_dimensions(self):
        return self._op.row_dimensions

    @property
    def column_dimensions(self):
        return self._op.column_dimensions
    
