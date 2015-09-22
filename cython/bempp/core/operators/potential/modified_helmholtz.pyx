__all__=['single_layer','double_layer']

from bempp.core.utils cimport shared_ptr
from bempp.core.utils cimport catch_exception
from bempp.core.utils cimport complex_double
from bempp.core.assembly.discrete_boundary_operator cimport c_DiscreteBoundaryOperator
from bempp.core.assembly.discrete_boundary_operator cimport RealDiscreteBoundaryOperator
from bempp.core.assembly.discrete_boundary_operator cimport ComplexDiscreteBoundaryOperator
from bempp.core.space cimport Space, c_Space
from bempp.core.utils cimport np_to_eigen_matrix_float64, Matrix
from bempp.core.utils cimport c_ParameterList, ParameterList
from cython.operator cimport dereference as deref

cimport numpy as _np
import numpy as _np

cdef extern from "bempp/core/operators/potential/py_potential_operators.hpp" namespace "Bempp":
    cdef shared_ptr[const c_DiscreteBoundaryOperator[complex_double]] modified_helmholtz_single_layer_potential_discrete_operator(
                const shared_ptr[const c_Space[double]]& space,
                const Matrix[double]& evaluation_points,
                complex_double wave_number,
                const c_ParameterList& parameterList) except +catch_exception

    cdef shared_ptr[const c_DiscreteBoundaryOperator[complex_double]] modified_helmholtz_double_layer_potential_discrete_operator(
                const shared_ptr[const c_Space[double]]& space,
                const Matrix[double]& evaluation_points,
                complex_double wave_number,
                const c_ParameterList& parameterList) except +catch_exception


def single_layer_ext(Space space not None,
        _np.ndarray evaluation_points not None, 
        double complex wave_number,
        ParameterList parameters not None):

        if not (evaluation_points.ndim==2 and evaluation_points.shape[0]==3):
            raise ValueError("Wrong format for input points")

        points = _np.require(evaluation_points,"double","F")

        cdef ComplexDiscreteBoundaryOperator op = ComplexDiscreteBoundaryOperator()

        op.impl_.assign(
             modified_helmholtz_single_layer_potential_discrete_operator(
                 space.impl_,np_to_eigen_matrix_float64(points),
                 complex_double(_np.real(wave_number), _np.imag(wave_number)),
                 deref(parameters.impl_)))
        return op

                
def double_layer_ext(Space space not None,
        _np.ndarray evaluation_points not None, 
        double complex wave_number,
        ParameterList parameters not None):

        if not (evaluation_points.ndim==2 and evaluation_points.shape[0]==3):
            raise ValueError("Wrong format for input points")

        points = _np.require(evaluation_points,"double","F")

        cdef ComplexDiscreteBoundaryOperator op = ComplexDiscreteBoundaryOperator()

        op.impl_.assign(
             modified_helmholtz_double_layer_potential_discrete_operator(
                 space.impl_,np_to_eigen_matrix_float64(points),
                 complex_double(_np.real(wave_number), _np.imag(wave_number)),
                 deref(parameters.impl_)))
        return op
