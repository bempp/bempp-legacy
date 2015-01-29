#cython: embedsignature=True

__all__=['single_layer','double_layer','adjoint_double_layer','hypersingular']

from bempp.utils cimport catch_exception, complex_float, complex_double
from bempp.assembly.boundary_operator cimport BoundaryOpVariants
from bempp.utils.parameter_list cimport c_ParameterList, ParameterList
from bempp.space.space cimport SpaceVariants,Space
from bempp.utils cimport shared_ptr
from libcpp.string cimport string
from numpy cimport dtype

cdef extern from "bempp/operators/py_operators.hpp" namespace "Bempp":
    BoundaryOpVariants c_modifiedHelmholtz3dSingleLayerBoundaryOperator[BASIS,RESULT](
            const c_ParameterList&,
            const SpaceVariants&,
            const SpaceVariants&,
            const SpaceVariants&,
            RESULT waveNumber,
            const string&,
            int) except+catch_exception


cdef extern from "bempp/operators/py_operators.hpp" namespace "Bempp":
    BoundaryOpVariants c_modifiedHelmholtz3dDoubleLayerBoundaryOperator[BASIS,RESULT](
            const c_ParameterList&,
            const SpaceVariants&,
            const SpaceVariants&,
            const SpaceVariants&,
            RESULT waveNumber,
            const string&,
            int) except+catch_exception

cdef extern from "bempp/operators/py_operators.hpp" namespace "Bempp":
    BoundaryOpVariants c_modifiedHelmholtz3dHypersingularBoundaryOperator[BASIS,RESULT](
            const c_ParameterList&,
            const SpaceVariants&,
            const SpaceVariants&,
            const SpaceVariants&,
            RESULT waveNumber,
            const string&,
            int) except+catch_exception

cdef extern from "bempp/operators/py_operators.hpp" namespace "Bempp":
    BoundaryOpVariants c_modifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator[BASIS,RESULT](
            const c_ParameterList&,
            const SpaceVariants&,
            const SpaceVariants&,
            const SpaceVariants&,
            RESULT waveNumber,
            const string&,
            int) except+catch_exception
