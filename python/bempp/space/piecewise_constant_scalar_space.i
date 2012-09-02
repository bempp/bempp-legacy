%{
#include "space/piecewise_constant_scalar_space.hpp"
%}

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
boost::shared_ptr<Space<BasisFunctionType> >
piecewiseConstantScalarSpace(Grid& grid)
{
    typedef PiecewiseConstantScalarSpace<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(grid));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(piecewiseConstantScalarSpace);
}

%pythoncode %{

def piecewiseConstantScalarSpace(context, grid):
    """Space of piecewise constant scalar functions"""
    name = 'PiecewiseConstantScalarSpace'
    return constructObjectTemplatedOnBasis(name, context.basisFunctionType(), grid)

%}

