__all__=['single_layer','double_layer']

from bempp.utils cimport shared_ptr
from bempp.assembly.discrete_boundary_operator cimport c_DiscreteBoundaryOperator
from bempp.space.space cimport SpaceVariants
from bempp.utils.armadillo cimport Mat
from bempp.utils.parameter_list cimport c_ParameterList, ParameterList
from bempp.utils.armadillo cimport Mat
from bempp.assembly.potential_operator cimport PotentialOperator
from bempp.assembly.discrete_boundary_operator cimport DiscreteBoundaryOperator
from bempp.common import global_parameters
from bempp.space.space cimport Space
from cython.operator cimport dereference as deref

import numpy as np
cimport numpy as np

cdef extern from "bempp/operators/py_operators.hpp" namespace "Bempp":
    cdef shared_ptr[const c_DiscreteBoundaryOperator[double]] py_laplace_single_layer_potential_discrete_operator(
                const SpaceVariants& space,
                const Mat[double]& evaluationPoints,
                const c_ParameterList& parameterList) 

    cdef shared_ptr[const c_DiscreteBoundaryOperator[double]] py_laplace_double_layer_potential_discrete_operator(
                const SpaceVariants& space,
                const Mat[double]& evaluationPoints,
                const c_ParameterList& parameterList) 


def single_layer(Space space,
        np.ndarray evaluation_points, ParameterList parameter_list=None):

        cdef int component_count = 1

        if parameter_list is None:
            parameter_list = global_parameters()

        if not (evaluation_points.ndim==2 and evaluation_points.shape[0]==3):
            raise ValueError("Wrong format for input points")

        points = np.require(evaluation_points,"double","F")
        cdef double[::1,:] point_data = points
        cdef shared_ptr[Mat[double]] arma_data = shared_ptr[Mat[double]](
                new Mat[double](&point_data[0,0],
                evaluation_points.shape[0],evaluation_points.shape[1],
                False,True))

        cdef DiscreteBoundaryOperator op = DiscreteBoundaryOperator()

        op._impl_float64_.assign(
             py_laplace_single_layer_potential_discrete_operator(
                 space.impl_,deref(arma_data),deref(parameter_list.impl_)))
        op._dtype = np.dtype('float64')

        return PotentialOperator(op,component_count,space,points)

                
def double_layer(Space space,
        np.ndarray evaluation_points, ParameterList parameter_list=None):

        cdef int component_count = 1

        if parameter_list is None:
            parameter_list = global_parameters()

        if not (evaluation_points.ndim==2 and evaluation_points.shape[0]==3):
            raise ValueError("Wrong format for input points")

        points = np.require(evaluation_points,"double","F")
        cdef double[::1,:] point_data = points
        cdef shared_ptr[Mat[double]] arma_data = shared_ptr[Mat[double]](
                new Mat[double](&point_data[0,0],
                evaluation_points.shape[0],evaluation_points.shape[1],
                False,True))

        cdef DiscreteBoundaryOperator op = DiscreteBoundaryOperator()

        op._impl_float64_.assign(
             py_laplace_double_layer_potential_discrete_operator(
                 space.impl_,deref(arma_data),deref(parameter_list.impl_)))
        op._dtype = np.dtype('float64')

        return PotentialOperator(op,component_count,space,points)

                     


