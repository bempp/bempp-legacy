%{
  #include "space/piecewise_constant_scalar_space.hpp"
%}

// TODO
// %include "space_docstrings.i"

namespace Bempp
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(PiecewiseConstantScalarSpace);
}

%include "space/piecewise_constant_scalar_space.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(PiecewiseConstantScalarSpace);
}

%pythoncode %{

def piecewiseConstantScalarSpace(grid, basisFunctionType='float64'):
    """Space of piecewise constant scalar functions"""
    name = 'PiecewiseConstantScalarSpace'
    return constructObjectTemplatedOnBasis(name, basisFunctionType, grid)

%}

