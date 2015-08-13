__all__=['single_layer','double_layer']

from bempp.utils cimport shared_ptr
from bempp.assembly.discrete_boundary_operator cimport c_DiscreteBoundaryOperator
from bempp.space.space cimport SpaceVariants
from bempp.utils.eigen cimport np_to_eigen_matrix_float64, Matrix
from bempp.utils.parameter_list cimport c_ParameterList, ParameterList
from bempp.assembly.potential_operator cimport PotentialOperator
from bempp.assembly.discrete_boundary_operator cimport DiscreteBoundaryOperator
from bempp import global_parameters
from bempp.space.space cimport Space
from cython.operator cimport dereference as deref

import numpy as _np
cimport numpy as _np

cdef extern from "bempp/operators/py_operators.hpp" namespace "Bempp":
    cdef shared_ptr[const c_DiscreteBoundaryOperator[double]] py_laplace_single_layer_potential_discrete_operator(
                const SpaceVariants& space,
                const Matrix[double]& evaluation_points,
                const c_ParameterList& parameterList) 

    cdef shared_ptr[const c_DiscreteBoundaryOperator[double]] py_laplace_double_layer_potential_discrete_operator(
                const SpaceVariants& space,
                const Matrix[double]& evaluation_points,
                const c_ParameterList& parameterList) 


def single_layer(Space space,
        _np.ndarray evaluation_points, ParameterList parameters=None):

        cdef int component_count = 1

        if parameters is None:
            parameters = global_parameters

        if not (evaluation_points.ndim==2 and evaluation_points.shape[0]==3):
            raise ValueError("Wrong format for input points")

        points = _np.require(evaluation_points,"double","F")

        cdef DiscreteBoundaryOperator op = DiscreteBoundaryOperator()

        op._impl_float64_.assign(
             py_laplace_single_layer_potential_discrete_operator(
                 space.impl_,np_to_eigen_matrix_float64(points),deref(parameters.impl_)))
        op._dtype = _np.dtype('float64')

        return PotentialOperator(op,component_count,space,points)

                
def double_layer(Space space,
        _np.ndarray evaluation_points, ParameterList parameters=None):

        cdef int component_count = 1

        if parameters is None:
            parameters = global_parameters

        if not (evaluation_points.ndim==2 and evaluation_points.shape[0]==3):
            raise ValueError("Wrong format for input points")

        points = _np.require(evaluation_points,"double","F")

        cdef DiscreteBoundaryOperator op = DiscreteBoundaryOperator()

        op._impl_float64_.assign(
             py_laplace_double_layer_potential_discrete_operator(
                 space.impl_,np_to_eigen_matrix_float64(points),deref(parameters.impl_)))
        op._dtype = _np.dtype('float64')

        return PotentialOperator(op,component_count,space,points)

                     


