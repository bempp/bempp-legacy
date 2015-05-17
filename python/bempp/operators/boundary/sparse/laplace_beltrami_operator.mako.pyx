#cython: embedsignature=True
<%
from data_types import dtypes, compatible_dtypes, ctypes
%>

from bempp.utils.parameter_list cimport c_ParameterList, ParameterList
from bempp.space.space cimport SpaceVariants,Space
from libcpp.string cimport string
from bempp.utils.enum_types cimport symmetry_mode
from bempp.assembly.boundary_operator cimport BoundaryOpVariants,GeneralBoundaryOperator
from cython.operator cimport dereference as deref
from bempp.utils.byte_conversion import convert_to_bytes
from bempp.utils.enum_types cimport symmetry_mode
from bempp.utils cimport complex_float, complex_double
from bempp import global_parameters

def laplace_beltrami(Space domain, Space range, Space dual_to_range,
        object label="", object symmetry="auto_symmetry", 
        parameters=None):

    cdef ParameterList local_parameters
    cdef GeneralBoundaryOperator bop 

    if not len({domain.dtype,range.dtype,dual_to_range.dtype})==1:
        raise ValueError("All spaces must have the same data type")


    if parameters is None:
        local_parameters = global_parameters
    else:
        if not isinstance(parameters,ParameterList):
            raise ValueError("parameters must be of type bempp.ParameterList")
        local_parameters = parameters


    basis_type = domain.dtype
    result_type = 'float64'

    bop = GeneralBoundaryOperator(basis_type,result_type,
            local_parameters,True)

    bop.impl_.assign(c_laplaceBeltramiOperator[double,double](
        deref(local_parameters.impl_),domain.impl_,range.impl_,
        dual_to_range.impl_,convert_to_bytes(label),
        symmetry_mode(convert_to_bytes(symmetry))))

    return bop





    










