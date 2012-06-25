%{
  #include "space/piecewise_linear_continuous_scalar_space.hpp"
%}

// TODO
// %include "space_docstrings.i"

namespace Bempp
{
BEMPP_FORWARD_DECLARE_CLASS_TEMPLATED_ON_BASIS(PiecewiseLinearContinuousScalarSpace);
}

%include "space/piecewise_linear_continuous_scalar_space.hpp"

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(PiecewiseLinearContinuousScalarSpace);
}

%pythoncode %{

def piecewiseLinearContinuousScalarSpace(grid, basisFunctionType='float64'):
    """Space of piecewise linear continuous scalar functions"""
    name = 'PiecewiseLinearContinuousScalarSpace'
    return constructObjectTemplatedOnBasis(name, basisFunctionType, grid)

%}

