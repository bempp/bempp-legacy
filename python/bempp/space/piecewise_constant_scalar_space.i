%{
  #include "space/piecewise_constant_scalar_space.hpp"
%}

// TODO
// %include "space_docstrings.i"



namespace Bempp
{

    BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_BASIS(PiecewiseConstantScalarSpace);
}

%include "space/piecewise_constant_scalar_space.hpp"


%pythoncode %{

def piecewiseConstantScalarSpace(grid,basisFunctionType='float64'):
    """Space of piecewise constant scalar functions"""
    name='PiecewiseConstantScalarSpace'
    return constructObjectTemplatedOnBasis(name, basisFunctionType,grid)
      %}
