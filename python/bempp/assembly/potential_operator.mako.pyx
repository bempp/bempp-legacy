from bempp.assembly.discrete_boundary_operator cimport DiscreteBoundaryOperator
from bempp.assembly.discrete_boundary_operator cimport DiscreteBoundaryOperatorBase
from bempp.space.space cimport Space
from bempp.assembly.grid_function cimport GridFunction
from bempp.assembly.grid_function import GridFunction
import numpy as np
cimport numpy as np



cdef class PotentialOperator:

    def __cinit__(self,DiscreteBoundaryOperatorBase op,
            int component_count,
            Space space, 
            np.ndarray evaluation_points):
        pass

    def __init__(self,DiscreteBoundaryOperatorBase op,
            int component_count,
            Space space, 
            np.ndarray evaluation_points):

        self._op = op
        self._component_count = component_count
        self._space = space
        self._evaluation_points = evaluation_points

    def evaluate(self, GridFunction g):

        res = self._op*g.coefficients

        return res.reshape(self._component_count,-1,order='F')

    def __is_compatible(self,PotentialOperator other):

        return (self.component_count==other.component_count and
                np.linalg.norm(self.evaluation_points-
                    other.evaluation_points,ord=np.inf)==0 and
                self.space.is_compatible(other.space))

    def __add__(self,PotentialOperator other):

        if not self.__is_compatible(other):
            raise ValueError("Potential Operator sum not possible")

        return PotentialOperator(
                self.discrete_operator+other.discrete_operator,
                self.component_count,
                self.space,self.evaluation_points)

    def __mul__(self,object obj):

        if not isinstance(self,PotentialOperator):
            return obj*self

        if np.isscalar(obj):
            return PotentialOperator(obj*self.discrete_operator,
                    self.component_count,
                    self.space,self.evaluation_points,)
        elif isinstance(obj,GridFunction):
            return self.evaluate(obj)
        else:
            raise NotImplementedError("Cannot multiply with object of type {0}".format(str(type(obj))))

    def __neg__(self):

        return self.__mul__(-1.0)

    def __sub__(self,PotentialOperator other):

        return self.__add__(-other)

    property space:

        def __get__(self):
            return self._space

    property component_count:

        def __get__(self):
            return self._component_count

    property evaluation_points:

        def __get__(self):
            return self._evaluation_points

    property discrete_operator:

        def __get__(self):
            return self._op

    property dtype:

        def __get__(self):
            return self.discrete_operator.dtype



