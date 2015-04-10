#cython: embedsignature=True

__all__=['single_layer','double_layer']

from bempp.utils cimport shared_ptr
from bempp.assembly.discrete_boundary_operator cimport c_DiscreteBoundaryOperator
from bempp.space.space cimport SpaceVariants
from bempp.utils.parameter_list cimport c_ParameterList, ParameterList
from bempp.utils.eigen cimport np_to_eigen_matrix_float64
from bempp.utils cimport complex_double
from bempp.assembly.potential_operator cimport PotentialOperator
from bempp.assembly.discrete_boundary_operator cimport DiscreteBoundaryOperator
from bempp.common import global_parameters
from bempp.space.space cimport Space
from cython.operator cimport dereference as deref

import numpy as np
cimport numpy as np

cdef extern from "bempp/operators/py_operators.hpp" namespace "Bempp":
    cdef shared_ptr[const c_DiscreteBoundaryOperator[complex_double]] py_helmholtz_single_layer_far_field_discrete_operator(
                const SpaceVariants& space,
                complex_double waveNumber,
                const c_ParameterList& parameterList) 

    cdef shared_ptr[const c_DiscreteBoundaryOperator[complex_double]] py_helmholtz_double_layer_far_field_discrete_operator(
                const SpaceVariants& space,
                complex_double waveNumber,
                const c_ParameterList& parameterList) 


def single_layer(Space space,
        np.ndarray evaluation_points, 
        object wave_number, ParameterList parameter_list=None):

        cdef int component_count = 1

        if parameter_list is None:
            parameter_list = global_parameters()

        if not (evaluation_points.ndim==2 and evaluation_points.shape[0]==3):
            raise ValueError("Wrong format for input points")

        points = np.require(evaluation_points,"double","F")

        cdef DiscreteBoundaryOperator op = DiscreteBoundaryOperator()
        cdef complex_double cpp_wave_number = complex_double(np.real(wave_number),
                np.imag(wave_number))

        op._impl_complex128_.assign(
             py_helmholtz_single_layer_far_field_discrete_operator(
                 space.impl_,deref(np_to_eigen_matrix_float64(points)),
                 cpp_wave_number,deref(parameter_list.impl_)))
        op._dtype = np.dtype('complex128')

        return PotentialOperator(op,component_count,space,points)

                
def double_layer(Space space,
        np.ndarray evaluation_points, 
        object wave_number, ParameterList parameter_list=None):

        cdef int component_count = 1

        if parameter_list is None:
            parameter_list = global_parameters()

        if not (evaluation_points.ndim==2 and evaluation_points.shape[0]==3):
            raise ValueError("Wrong format for input points")

        points = np.require(evaluation_points,"double","F")

        cdef DiscreteBoundaryOperator op = DiscreteBoundaryOperator()
        cdef complex_double cpp_wave_number = complex_double(np.real(wave_number),
                np.imag(wave_number))

        op._impl_complex128_.assign(
             py_helmholtz_double_layer_far_field_discrete_operator(
                 space.impl_,deref(np_to_eigen_matrix_float64(points)),
                 cpp_wave_number,deref(parameter_list.impl_)))
        op._dtype = np.dtype('complex128')

        return PotentialOperator(op,component_count,space,points)

                     


