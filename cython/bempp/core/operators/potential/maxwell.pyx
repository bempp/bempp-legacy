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

cdef extern from "bempp/operators/maxwell_operators.hpp" namespace "Bempp":
    cdef shared_ptr[const c_DiscreteBoundaryOperator[complex_double]] maxwell_electric_field_potential_operator "Bempp::electricFieldPotentialOperator<double>"(
                const shared_ptr[const c_Space[double]]& space,
                const Matrix[double]& evaluation_points,
                complex_double wave_number,
                const c_ParameterList& parameterList) except +catch_exception

    cdef shared_ptr[const c_DiscreteBoundaryOperator[complex_double]] maxwell_magnetic_field_potential_operator "Bempp::magneticFieldPotentialOperator"(
                const shared_ptr[const c_Space[double]]& space,
                const Matrix[double]& evaluation_points,
                complex_double wave_number,
                const c_ParameterList& parameterList) except +catch_exception


def electric_field_ext(Space space not None,
        _np.ndarray evaluation_points not None, 
        double complex wave_number,
        ParameterList parameters not None):

        if not (evaluation_points.ndim==2 and evaluation_points.shape[0]==3):
            raise ValueError("Wrong format for input points")

        points = _np.require(evaluation_points,"double","F")

        cdef ComplexDiscreteBoundaryOperator op = ComplexDiscreteBoundaryOperator()

        op.impl_.assign(
             maxwell_electric_field_potential_operator(
                 space.impl_,np_to_eigen_matrix_float64(points),
                 complex_double(_np.real(wave_number), _np.imag(wave_number)),
                 deref(parameters.impl_)))
        return op

                
def magnetic_field_ext(Space space not None,
        _np.ndarray evaluation_points not None, 
        double complex wave_number,
        ParameterList parameters not None):

        if not (evaluation_points.ndim==2 and evaluation_points.shape[0]==3):
            raise ValueError("Wrong format for input points")

        points = _np.require(evaluation_points,"double","F")

        cdef ComplexDiscreteBoundaryOperator op = ComplexDiscreteBoundaryOperator()

        op.impl_.assign(
             maxwell_magnetic_field_potential_operator(
                 space.impl_,np_to_eigen_matrix_float64(points),
                 complex_double(_np.real(wave_number), _np.imag(wave_number)),
                 deref(parameters.impl_)))
        return op
