%{
  #include "space/piecewise_constant_scalar_space.hpp"
%}

// TODO
// %include "space_docstrings.i"

%shared_ptr(Bempp::PiecewiseConstantScalarSpace<float>);
%shared_ptr(Bempp::PiecewiseConstantScalarSpace<double>);
%shared_ptr(Bempp::PiecewiseConstantScalarSpace<std::complex<float> >);
%shared_ptr(Bempp::PiecewiseConstantScalarSpace<std::complex<double> >);

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

def piecewiseConstantScalarSpace(context, grid):
    """Space of piecewise constant scalar functions"""
    name = 'PiecewiseConstantScalarSpace'
    return constructObjectTemplatedOnBasis(name, context.basisFunctionType(), grid)

%}

