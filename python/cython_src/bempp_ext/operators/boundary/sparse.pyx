from bempp_ext.assembly.abstract_boundary_operator cimport ElementaryLocalOperator
from bempp_ext.assembly.abstract_boundary_operator cimport c_ElementaryLocalOperator
from bempp_ext.space.space cimport c_Space, Space
from bempp_ext.utils cimport shared_ptr
from bempp_ext.utils.enum_types cimport SymmetryMode, symmetry_mode
from bempp_ext.utils cimport ParameterList, c_ParameterList
from libcpp.string cimport string
from cython.operator cimport dereference as deref


cdef extern from "bempp_ext/operators/boundary/py_boundary_operators.hpp" namespace "Bempp":
    shared_ptr[const c_ElementaryLocalOperator] identity_operator(
            const c_ParameterList&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            string label, SymmetryMode symmetry)
    shared_ptr[const c_ElementaryLocalOperator] maxwell_identity_operator(
            const c_ParameterList&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            string label, SymmetryMode symmetry)
    shared_ptr[const c_ElementaryLocalOperator] laplace_beltrami_operator(
            const c_ParameterList&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            string label, SymmetryMode symmetry)
    shared_ptr[const c_ElementaryLocalOperator] curl_value_local_operator(
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            shared_ptr[const c_Space[double]]&,
            int component)


def _convert_to_bytes(s):
    res = s
    try:
        if not isinstance(s,bytes):
            res = res.encode('UTF-8')
    except:
        raise ValueError('String type expected.')
    return res


def identity_ext(
        ParameterList parameters,
        Space domain,
        Space range,
        Space dual_to_range,
        object label='', object symmetry='no_symmetry'):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(identity_operator(
        deref(parameters.impl_),domain.impl_, range.impl_, dual_to_range.impl_,
        _convert_to_bytes(label), symmetry_mode(symmetry)))
    return op

def maxwell_identity_ext(
        ParameterList parameters,
        Space domain,
        Space range,
        Space dual_to_range,
        object label='', object symmetry='no_symmetry'):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(maxwell_identity_operator(
        deref(parameters.impl_),domain.impl_, range.impl_, dual_to_range.impl_,
        _convert_to_bytes(label), symmetry_mode(symmetry)))
    return op

def laplace_beltrami_ext(
        ParameterList parameters,
        Space domain,
        Space range,
        Space dual_to_range,
        object label='', object symmetry='no_symmetry'):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(laplace_beltrami_operator(
        deref(parameters.impl_),domain.impl_, range.impl_, dual_to_range.impl_,
        _convert_to_bytes(label), symmetry_mode(symmetry)))
    return op

def curl_value_ext(Space domain, Space range, Space dual_to_range,
                                  int component):

    cdef ElementaryLocalOperator op = ElementaryLocalOperator()
    op.impl_.assign(curl_value_local_operator(domain.impl_, range.impl_, dual_to_range.impl_,
                    component))
    return op