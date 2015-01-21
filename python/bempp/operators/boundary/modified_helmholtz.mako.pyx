#cython: embedsignature=True
<%
from data_types import dtypes, compatible_dtypes, ctypes, scalar_cython_type
op_names = [('single_layer','c_modifiedHelmholtz3dSingleLayerBoundaryOperator',
             'Return the modified Helmholtz single layer boundary operator.'),
            ('double_layer','c_modifiedHelmholtz3dDoubleLayerBoundaryOperator',
             'Return the modified Helmholtz double layer boundary operator.'),
            ('adjoint_double_layer','c_modifiedHelmholtz3dAdjointDoubleLayerBoundaryOperator',
             'Return the modified Helmholtz adjoint double layer boundary operator.'),
            ('hypersingular','c_modifiedHelmholtz3dHypersingularBoundaryOperator',
              'Return the modified Helmholtz hypersingular boundary operator.')]

%>

from bempp.utils.parameter_list cimport c_ParameterList, ParameterList
from bempp.space.space cimport SpaceVariants,Space
from libcpp.string cimport string
from bempp.utils.enum_types cimport symmetry_mode
from bempp.assembly.boundary_operator cimport GeneralBoundaryOperator,DenseBoundaryOperator,BoundaryOpVariants
from cython.operator cimport dereference as deref
from bempp.utils.byte_conversion import convert_to_bytes
from bempp.utils.enum_types cimport symmetry_mode
from bempp.utils cimport complex_float, complex_double
from bempp.common import global_parameters
import numpy as np
cimport numpy as np

% for pyname,c_name,help_text in op_names:
def ${pyname}(Space domain, Space range, Space dual_to_range,
        object wave_number,
        object label="", object symmetry="auto_symmetry", 
        object result_type=None,
        object parameter_list=None):
    """

    ${help_text}

    """

    cdef ParameterList parameters
    cdef GeneralBoundaryOperator bop 


% for pyresult,cyresult in dtypes.items():
    cdef ${scalar_cython_type(cyresult)} cy_wave_number_${pyresult}
    cdef ${cyresult} c_wave_number_${pyresult}
% endfor

    if not len({domain.dtype,range.dtype,dual_to_range.dtype})==1:
        raise ValueError("All spaces must have the same data type")

    if parameter_list is None:
        parameters = global_parameters()
    else:
        if not isinstance(parameter_list,ParameterList):
            raise ValueError("parameter_list must be of type bempp.ParameterList")
        parameters = parameter_list

    actual_result_type = None
    basis_type = domain.dtype

    if result_type is None:
        if np.iscomplexobj(wave_number):
            if basis_type=='float32':
                actual_result_type = 'complex64'
            elif basis_type=='float64':
                actual_result_type = 'complex128'
            else:
                raise ValueError("Only 'float32' and 'float64' are supported as basis types.")
        else:
            if basis_type=='float32':
                actual_result_type = 'float32'
            elif basis_type=='float64':
                actual_result_type = 'float64'
            else:
                raise ValueError("Only 'float32' and 'float64' are supported as basis types.")
    else:
        actual_result_type = result_type



    bop = GeneralBoundaryOperator(basis_type=basis_type,result_type=actual_result_type)


% for pybasis, cybasis in dtypes.items():
%     for pyresult, cyresult in dtypes.items():
%         if pyresult in compatible_dtypes[pybasis]:


    if basis_type=="${pybasis}" and actual_result_type=="${pyresult}":

        cy_wave_number_${pyresult} = wave_number
        c_wave_number_${pyresult} = deref(<${cyresult}*>&cy_wave_number_${pyresult})

        bop.impl_.assign(${c_name}[${cybasis},${cyresult}](
            deref(parameters.impl_),domain.impl_,range.impl_,
            dual_to_range.impl_, c_wave_number_${pyresult}, convert_to_bytes(label),
            symmetry_mode(convert_to_bytes(symmetry))))
        if 'boundaryOperatorAssemblyType' in parameters and parameters['boundaryOperatorAssemblyType']=='dense': 
            return DenseBoundaryOperator(bop)
%           endif
%       endfor
% endfor

    raise ValueError("Wrong basis_type or result_type")


% endfor
