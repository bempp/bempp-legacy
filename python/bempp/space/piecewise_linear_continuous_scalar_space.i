%{
  #include "space/piecewise_linear_continuous_scalar_space.hpp"
%}

// TODO
// %include "space_docstrings.i"

%shared_ptr(Bempp::PiecewiseLinearContinuousScalarSpace<float>);
%shared_ptr(Bempp::PiecewiseLinearContinuousScalarSpace<double>);
%shared_ptr(Bempp::PiecewiseLinearContinuousScalarSpace<std::complex<float> >);
%shared_ptr(Bempp::PiecewiseLinearContinuousScalarSpace<std::complex<double> >);

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

def piecewiseLinearContinuousScalarSpace(context, grid):
    """Space of piecewise linear continuous scalar functions"""
    name = 'PiecewiseLinearContinuousScalarSpace'
    return constructObjectTemplatedOnBasis(name, context.basisFunctionType(), grid)

%}

