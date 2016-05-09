import numpy as _np
from bempp.api.utils.linear_operator import LinearOperator


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


class BlockedOperator(object):

    def __init__(self, m, n):

        self._m = m
        self._n = n

        self._operators = _np.empty((m, n), dtype=_np.object)
        self._rows = m * [False]
        self._cols = n * [False]

    def __getitem__(self, key):

        return self._operators[key]

    def __setitem__(self, key, operator):

        self._operators[key] = operator
        self._rows[key[0]] = True
        self._cols[key[1]] = True

    def _fill_complete(self):

        return (False not in self._cols) and (False not in self._rows)

    def weak_form(self):

        if not self._fill_complete():
            raise ValueError(
                "Each row and column must have at least one operator")

        discrete_operator = BlockedDiscreteOperator(self._m, self._n)

        for i in range(self._m):
            for j in range(self._n):
                if self._operators[i, j] is not None:
                    discrete_operator[i, j] = self._operators[i, j].weak_form()

        return discrete_operator

    def strong_form(self, mode='simple'):

        if not self._fill_complete():
            raise ValueError(
                "Each row and column must have at least one operator")

        if mode == 'full':
            discrete_operator = BlockedDiscreteOperator(self._m, self._n)

            for i in range(self._m):
                for j in range(self._n):
                    if self._operators[i, j] is not None:
                        discrete_operator[i, j] = self._operators[
                            i, j].strong_form()
            return discrete_operator

        elif mode == 'simple':
            from bempp.api import InverseSparseDiscreteBoundaryOperator
            from bempp.api.operators.boundary.sparse import identity

            blocked = BlockedDiscreteOperator(self.ndims[0], self.ndims[0])

            for i in range(self.ndims[0]):
                op = None
                for j in range(self.ndims[1]):
                    if self[i, j] is not None:
                        op = self[i, j]
                        break

                blocked[i, i] = InverseSparseDiscreteBoundaryOperator(
                    identity(op.range, op.dual_to_range,
                             op.dual_to_range).weak_form())
            return blocked * self.weak_form()
        else:
            raise ValueError(
                "Unknown value for 'mode'. Allowed values are 'simple' and 'full'")

    def __add__(self, other):

        if self.ndims != other.ndims:
            return ValueError("Both blocked operators must have the same dimensions.")

        blocked_operator = BlockedOperator(self.ndims[0], self.ndims[1])

        for i in range(self.ndims[0]):
            for j in range(self.ndims[1]):
                if other[i, j] is None:
                    blocked_operator[i, j] = self[i, j]
                elif self[i, j] is None:
                    blocked_operator[i, j] = other[i, j]
                else:
                    blocked_operator[i, j] = self[i, j] + other[i, j]
        return blocked_operator

    def __neg__(self):
        return self.__mul__(-1)

    def __sub__(self, other):
        return self.__add__(-other)

    def __mul__(self, other):

        import numpy as np
        if np.isscalar(other):

            blocked_operator = BlockedOperator(self.ndims[0], self.ndims[1])
            for i in range(self.ndims[0]):
                for j in range(self.ndims[1]):
                    if self[i, j] is not None:
                        blocked_operator[i, j] = other * self[i, j]
            return blocked_operator
        elif isinstance(other, BlockedOperator):
            if self.ndims[1] != other.ndims[0]:
                return ValueError("Dimensions are not compatible.")
            blocked_operator = BlockedOperator(self.ndims[0], other.ndims[1])
            for i in range(self.ndims[0]):
                for j in range(other.ndims[1]):
                    for k in range(self.ndims[1]):
                        blocked_operator[i, j] = _sum(blocked_operator[i, j],
                                                      _prod(self[i, k], other[k, j]))
            return blocked_operator
        else:
            return NotImplementedError

    def __rmul__(self, other):

        import numpy as np
        if np.isscalar(other):
            return self.__mul__(other)
        else:
            return NotImplementedError

    def get_ndims(self):

        return (self._m, self._n)

    ndims = property(get_ndims)


class BlockedDiscreteOperator(LinearOperator):

    def __init__(self, m, n):

        self._m = m
        self._n = n
        self._operators = _np.empty((m, n), dtype=_np.object)
        self._rows = _np.zeros(m, dtype=int)
        self._cols = _np.zeros(n, dtype=int)

    def __getitem__(self, key):

        return self._operators[key]

    def __setitem__(self, key, operator):

        if self._rows[key[0]] != 0:
            if operator.shape[0] != self._rows[key[0]]:
                raise ValueError("Incompatible number of rows")
        else:
            self._rows[key[0]] = operator.shape[0]

        if self._cols[key[1]] != 0:
            if operator.shape[1] != self._cols[key[1]]:
                raise ValueError("Incompatible number of columns")
        else:
            self._cols[key[1]] = operator.shape[1]
        self._operators[key] = operator

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

        for i in range(self._m):
            col_dim = 0
            local_res = res[row_dim:row_dim + self._rows[i], :]
            for j in range(self._n):
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

    def _get_shape(self):
        return (_np.sum(self._rows), _np.sum(self._cols))

    def _get_dtype(self):

        from bempp.api.utils.data_types import combined_type

        d = 'float64'
        for obj in self._operators.ravel():
            if obj is not None:
                d = combined_type(d, obj.dtype)

        return d

    def _get_ndims(self):

        return (self._m, self._n)

    def _get_row_dimensions(self):
        return self._rows

    def _get_column_dimensions(self):
        return self._cols

    def _as_matrix(self):
        if not self._fill_complete():
            raise ValueError("Not all rows or columns contain operators.")
        rows = []
        for i in range(self._m):
            row = []
            for j in range(self._n):
                if self[i, j] is None:
                    row.append(_np.zeros((self._rows[i], self._cols[j])))
                else:
                    row.append(self[i, j].as_matrix())
            rows.append(_np.hstack(row))
        return _np.vstack(rows)

    shape = property(_get_shape)
    dtype = property(_get_dtype)
    ndims = property(_get_ndims)
    row_dimensions = property(_get_row_dimensions)
    column_dimensions = property(_get_column_dimensions)
