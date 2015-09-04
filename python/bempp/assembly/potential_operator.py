"""Definition of potential operators."""


class PotentialOperator:
    def __init__(self, op, component_count, space, evaluation_points):

        self._op = op
        self._component_count = component_count
        self._space = space
        self._evaluation_points = evaluation_points

    def evaluate(self, grid_fun):

        res = self._op * grid_fun.coefficients
        return res.reshape(self._component_count, -1, order='F')

    def __is_compatible(self, other):
        import numpy as np

        return (self.component_count == other.component_count and
                np.linalg.norm(self.evaluation_points -
                               other.evaluation_points, ord=np.inf) == 0 and
                self.space.is_compatible(other.space))

    def __add__(self, other):

        if not self.__is_compatible(other):
            raise ValueError("Potential operators not compatible.")

        return PotentialOperator(
            self.discrete_operator + other.discrete_operator,
            self.component_count,
            self.space, self.evaluation_points)

    def __mul__(self, obj):
        import numpy as np
        from bempp import GridFunction

        if not isinstance(self, PotentialOperator):
            return obj * self

        if np.isscalar(obj):
            return PotentialOperator(obj * self.discrete_operator,
                                     self.component_count,
                                     self.space, self.evaluation_points, )
        elif isinstance(obj, GridFunction):
            return self.evaluate(obj)
        else:
            raise NotImplementedError("Cannot multiply with object of type {0}".format(str(type(obj))))

    def __neg__(self):

        return self.__mul__(-1.0)

    def __sub__(self, other):

        return self.__add__(-other)

    @property
    def space(self):
        return self._space

    @property
    def component_count(self):
        return self._component_count

    @property
    def evaluation_points(self):
        return self._evaluation_points

    @property
    def discrete_operator(self):
        return self._op

    @property
    def dtype(self):
        return self._op.dtype
