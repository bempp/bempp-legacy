%{
  #include "space/piecewise_linear_continuous_scalar_space.hpp"
%}

// TODO
// %include "space_docstrings.i"


namespace Bempp
{

    BEMPP_PYTHON_DECLARE_CLASS_TEMPLATED_ON_BASIS(PiecewiseConstantScalarSpace);
}

%include "space/piecewise_linear_continuous_scalar_space.hpp"





%pythoncode %{

def piecewiseLinearContinuousScalarSpace(grid,basisFunctionType='float64'):
    """Space of piecewise linear continuous scalar functions"""
    name='PiecewiseLinearContinuousScalarSpace'
    return constructObjectTemplatedOnBasis(name, basisFunctionType,grid)

	  %}

