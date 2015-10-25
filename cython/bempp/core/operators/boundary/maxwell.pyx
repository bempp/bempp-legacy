from bempp.core.assembly.abstract_boundary_operator cimport RealElementaryIntegralOperator
from bempp.core.assembly.abstract_boundary_operator cimport ComplexElementaryIntegralOperator
from bempp.core.assembly.abstract_boundary_operator cimport c_RealElementaryIntegralOperator
from bempp.core.assembly.abstract_boundary_operator cimport c_ComplexElementaryIntegralOperator
from bempp.core.space.space cimport c_Space, Space
from bempp.core.utils cimport shared_ptr
from bempp.core.utils.enum_types cimport SymmetryMode, symmetry_mode
from bempp.core.utils cimport ParameterList, c_ParameterList
from bempp.core.utils cimport complex_double
from libcpp.string cimport string
from cython.operator cimport dereference as deref

import numpy as np
cimport numpy as np

cdef extern from "bempp/operators/maxwell_operators.hpp" namespace "Bempp":
    shared_ptr[const c_ComplexElementaryIntegralOperator] maxwell_electric_field "Bempp::maxwellElectricFieldBoundaryOperator<double>"(
            const c_ParameterList&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            complex_double,
            string label, SymmetryMode symmetry)
    shared_ptr[const c_ComplexElementaryIntegralOperator] maxwell_magnetic_field "Bempp::maxwellMagneticFieldBoundaryOperator<double>"(
            const c_ParameterList&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            complex_double,
            string label, SymmetryMode symmetry)

def _convert_to_bytes(s):
    res = s
    try:
        if not isinstance(s,bytes):
            res = res.encode('UTF-8')
    except:
        raise ValueError('String type expected.')
    return res


def electric_field_ext(
        ParameterList parameters,
        Space domain,
        Space range,
        Space dual_to_range,
        double complex wave_number,
        object label='', object symmetry='no_symmetry'):

    cdef ComplexElementaryIntegralOperator op = ComplexElementaryIntegralOperator()
    op.impl_.assign(maxwell_electric_field(
        deref(parameters.impl_),domain.impl_, range.impl_, dual_to_range.impl_,
        complex_double(np.real(wave_number),np.imag(wave_number)),
        _convert_to_bytes(label), symmetry_mode(_convert_to_bytes(symmetry))))
    return op

def magnetic_field_ext(
        ParameterList parameters,
        Space domain,
        Space range,
        Space dual_to_range,
        double complex wave_number,
        object label='', object symmetry='no_symmetry'):

    cdef ComplexElementaryIntegralOperator op = ComplexElementaryIntegralOperator()
    op.impl_.assign(maxwell_magnetic_field(
        deref(parameters.impl_),domain.impl_, range.impl_, dual_to_range.impl_,
        complex_double(np.real(wave_number),np.imag(wave_number)),
        _convert_to_bytes(label), symmetry_mode(_convert_to_bytes(symmetry))))
    return op
