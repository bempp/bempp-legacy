%{
  #include "space/piecewise_constant_scalar_space.hpp"
%}

// TODO
// %include "space_docstrings.i"

namespace Bempp
{
BEMPP_PYTHON_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(PiecewiseConstantScalarSpace);
BEMPP_PYTHON_EXTEND_CONCRETE_CLASS_TEMPLATED_ON_BASIS(PiecewiseConstantScalarSpace);
}

%include "space/piecewise_constant_scalar_space.hpp"

namespace Bempp
{
BEMPP_PYTHON_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(PiecewiseConstantScalarSpace);
}

%pythoncode %{

def piecewiseConstantScalarSpace(grid,basisFunctionType='float64'):
    """Space of piecewise constant scalar functions"""
    name='PiecewiseConstantScalarSpace'
    return constructObjectTemplatedOnBasis(name, basisFunctionType,grid)

%}

