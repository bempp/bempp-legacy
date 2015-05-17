#cython: embedsignature=True

__all__=['single_layer','double_layer','adjoint_double_layer','hypersingular']

from bempp.utils.parameter_list cimport ParameterList
from bempp.space.space cimport Space
from bempp.operators.boundary import modified_helmholtz as _modified_helmholtz

<% ops = [('single_layer','Return the Helmholtz single layer boundary operator.'),
          ('double_layer','Return the Helmholtz double layer boundary operator.'),
          ('adjoint_double_layer','Return the Helmholtz adjoint double layer boundary operator.'),
          ('hypersingular','Return the Helmholtz hypersingular boundary operator.')]
    
%>

% for op,help_text in ops:
def ${op}(Space domain, Space range, Space dual_to_range,
        object wave_number,
        object label="",
        object symmetry="auto_symmetry",
        object parameters=None):
    """ 

    ${help_text}

    """

    return _modified_helmholtz.${op}(domain,range,dual_to_range,
            wave_number/1j,label,symmetry,parameters)


% endfor




