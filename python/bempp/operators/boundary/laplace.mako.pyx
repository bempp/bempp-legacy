#cython: embedsignature=True
<%
from data_types import dtypes, compatible_dtypes, ctypes

op_names = [('single_layer','c_laplace3dSingleLayerBoundaryOperator',
             'Return the Laplace single layer boundary operator.'),
            ('double_layer','c_laplace3dDoubleLayerBoundaryOperator',
             'Return the Laplace double layer boundary operator.'),
            ('adjoint_double_layer','c_laplace3dAdjointDoubleLayerBoundaryOperator',
             'Return the Laplace adjoint double layer boundary operator.'),
            ('hypersingular','c_laplace3dHypersingularBoundaryOperator',
              'Return the Laplace hypersingular boundary operator.')]

%>

__all__=['single_layer','double_layer','adjoint_double_layer','hypersingular']

from bempp.utils.parameter_list cimport c_ParameterList, ParameterList
from bempp.space.space cimport SpaceVariants,Space
from libcpp.string cimport string
from bempp.utils.enum_types cimport symmetry_mode
from bempp.assembly.boundary_operator cimport GeneralBoundaryOperator,BoundaryOpVariants
from bempp.assembly.boundary_operator cimport DenseBoundaryOperator
from cython.operator cimport dereference as deref
from bempp.utils.byte_conversion import convert_to_bytes
from bempp.utils.enum_types cimport symmetry_mode
from bempp.utils cimport complex_float, complex_double
from bempp.common import global_parameters

% for pyname,c_name,help_text in op_names:
def ${pyname}(Space domain, Space range, Space dual_to_range,
        object label="", object symmetry="auto_symmetry", 
        object parameter_list=None):
    """

    ${help_text}

    """

    cdef ParameterList parameters
    cdef GeneralBoundaryOperator bop 

    if not len({domain.dtype,range.dtype,dual_to_range.dtype})==1:
        raise ValueError("All spaces must have the same data type")

    if parameter_list is None:
        parameters = global_parameters()
    else:
        if not isinstance(parameter_list,ParameterList):
            raise ValueError("parameter_list must be of type bempp.ParameterList")
        parameters = parameter_list

    result_type = 'float64'

    basis_type = domain.dtype
    bop = GeneralBoundaryOperator(basis_type,result_type,
            parameters)

    bop.impl_.assign(
            ${c_name}[double,double](
            deref(parameters.impl_),domain.impl_,range.impl_,
            dual_to_range.impl_,convert_to_bytes(label),
            symmetry_mode(convert_to_bytes(symmetry))))

    return bop

%endfor
