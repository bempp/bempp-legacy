from bempp_ext.assembly.abstract_boundary_operator cimport RealElementaryIntegralOperator
from bempp_ext.assembly.abstract_boundary_operator cimport ComplexElementaryIntegralOperator
from bempp_ext.assembly.abstract_boundary_operator cimport c_RealElementaryIntegralOperator
from bempp_ext.assembly.abstract_boundary_operator cimport c_ComplexElementaryIntegralOperator
from bempp_ext.space.space cimport c_Space, Space
from bempp_ext.utils cimport shared_ptr
from bempp_ext.utils.enum_types cimport SymmetryMode
from bempp_ext.utils cimport ParameterList, c_ParameterList
from libcpp.string cimport string
from cython.operator cimport dereference as deref


cdef extern from "bempp_ext/operators/boundary/py_boundary_operators.hpp" namespace "Bempp":
    shared_ptr[const c_RealElementaryIntegralOperator] laplace_single_layer(
            const c_ParameterList&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            string label, SymmetryMode symmetry)
    shared_ptr[const c_RealElementaryIntegralOperator] laplace_double_layer(
            const c_ParameterList&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            string label, SymmetryMode symmetry)
    shared_ptr[const c_RealElementaryIntegralOperator] laplace_adjoint_double_layer(
            const c_ParameterList&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            string label, SymmetryMode symmetry)
    shared_ptr[const c_RealElementaryIntegralOperator] laplace_hypersingular(
            const c_ParameterList&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            string label, SymmetryMode symmetry)

def _convert_to_bytes(s):
    res = s
    try:
        if not isinstance(s,bytes):
            res = res.encode('UTF-8')
    except:
        raise ValueError('String type expected.')
    return res


def single_layer_ext(
        ParameterList parameters,
        Space domain,
        Space range,
        Space dual_to_range,
        object label='', object symmetry='no_symmetry'):

    cdef RealElementaryIntegralOperator op = RealElementaryIntegralOperator()
    op.impl_.assign(laplace_single_layer(
        deref(parameters.impl_),domain.impl_, range.impl_, dual_to_range.impl_,
        _convert_to_bytes(label), symmetry))
    return op

def double_layer_ext(
        ParameterList parameters,
        Space domain,
        Space range,
        Space dual_to_range,
        object label='', object symmetry='no_symmetry'):

    cdef RealElementaryIntegralOperator op = RealElementaryIntegralOperator()
    op.impl_.assign(laplace_double_layer(
        deref(parameters.impl_),domain.impl_, range.impl_, dual_to_range.impl_,
        _convert_to_bytes(label), symmetry))
    return op
        
def adjoint_double_layer_ext(
        ParameterList parameters,
        Space domain,
        Space range,
        Space dual_to_range,
        object label='', object symmetry='no_symmetry'):

    cdef RealElementaryIntegralOperator op = RealElementaryIntegralOperator()
    op.impl_.assign(laplace_adjoint_double_layer(
        deref(parameters.impl_),domain.impl_, range.impl_, dual_to_range.impl_,
        _convert_to_bytes(label), symmetry))
    return op
    
def hypersingular_ext(
        ParameterList parameters,
        Space domain,
        Space range,
        Space dual_to_range,
        object label='', object symmetry='no_symmetry'):

    cdef RealElementaryIntegralOperator op = RealElementaryIntegralOperator()
    op.impl_.assign(laplace_hypersingular(
        deref(parameters.impl_),domain.impl_, range.impl_, dual_to_range.impl_,
        _convert_to_bytes(label), symmetry))
    return op

