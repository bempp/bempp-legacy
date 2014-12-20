from bempp.utils.parameter_list cimport c_ParameterList, ParameterList
from bempp.space.space cimport SpaceVariants,Space
from libcpp.string cimport string
from bempp.utils.enum_types cimport symmetry_mode
from bempp.assembly.boundary_operator cimport BoundaryOperator,BoundaryOpVariants
from cython.operator cimport dereference as deref
from bempp.operators.sparse.identity_operator cimport c_identityOperator
from bempp.utils.byte_conversion import convert_to_bytes
from bempp.utils.enum_types cimport symmetry_mode

def identity_operator(ParameterList parameterList,
        Space domain, Space range, Space dual_to_range,
        object label="", object symmetry="auto_symmetry", 
        object basis_type="float64", object result_type="float64"):


    cdef BoundaryOperator bop = BoundaryOperator(basis_type='float64',result_type='float64')
    bop.impl_.assign(c_identityOperator[double,double](deref(parameterList.impl_),
            domain.impl_,range.impl_,dual_to_range.impl_,
            convert_to_bytes(label),symmetry_mode(convert_to_bytes(symmetry))))
    return bop





    










