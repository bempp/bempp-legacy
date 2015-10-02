import numpy as _np
from bempp.api.utils.linear_operator import LinearOperator


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
            raise ValueError("Each row and column must have at least one operator")

        discrete_operator = BlockedDiscreteOperator(self._m, self._n)

        for i in range(self._m):
            for j in range(self._n):
                if self._operators[i, j] is not None:
                    discrete_operator[i, j] = self._operators[i, j].weak_form()

        return discrete_operator

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

        super(BlockedDiscreteOperator, self).__init__(shape=self._get_shape(), dtype=self._get_dtype())

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

    def matvec(self, x):

        from .data_types import combined_type

        if x.ndim == 1:
            x_new = _np.expand_dims(x, 1)
            return self.matvec(x_new).ravel()

        if not self._fill_complete():
            raise ValueError("Not all rows or columns contain operators.")

        row_dim = 0
        res = _np.zeros((self.shape[0], x.shape[1]), dtype=combined_type(self.dtype, x.dtype))

        for i in range(self._m):
            col_dim = 0
            local_res = res[row_dim:row_dim + self._rows[i], :]
            for j in range(self._n):
                local_x = x[col_dim:col_dim + self._cols[j], :]
                if self._operators[i, j] is not None:
                    op_is_complex = _np.iscomplexobj(self._operators[i, j].dtype.type(1))
                    if _np.iscomplexobj(x) and not op_is_complex:
                        local_res[:] += (self._operators[i, j] * _np.real(local_x) +
                                         1j * self._operators[i, j] * _np.imag(local_x))
                    else:
                        local_res[:] += self._operators[i, j].dot(local_x)
                col_dim += self._cols[j]
            row_dim += self._rows[i]
        return res

    def matmat(self, x):

        return self.matvec(x)

    def __mul__(self, x):

        if _np.isscalar(x):
            res = BlockedDiscreteOperator(self.ndims[0], self.ndims[1])
            for i in range(self.ndims[0]):
                for j in range(self.ndims[1]):
                    op = self[i, j]
                    if op is not None: res[i, j] = x * op
            return res

        if isinstance(x, _np.ndarray):
            return self.matvec(x)

        raise NotImplementedError("Cannot multiply with object of type {0}".format(str(type(x))))

    def __rmul__(self, x):

        return self * x

    def __neg__(self):

        return self.__mul__(-1)

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

    ndims = property(_get_ndims)
    row_dimensions = property(_get_row_dimensions)
    column_dimensions = property(_get_column_dimensions)
