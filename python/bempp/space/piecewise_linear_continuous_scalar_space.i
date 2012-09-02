%{
#include "space/piecewise_linear_continuous_scalar_space.hpp"
%}

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
boost::shared_ptr<Space<BasisFunctionType> >
piecewiseLinearContinuousScalarSpace(Grid& grid)
{
    typedef PiecewiseLinearContinuousScalarSpace<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(grid));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(piecewiseLinearContinuousScalarSpace);
}

%pythoncode %{

def piecewiseLinearContinuousScalarSpace(context, grid):
    """Space of piecewise linear continuous scalar functions"""
    name = 'PiecewiseLinearContinuousScalarSpace'
    return constructObjectTemplatedOnBasis(name, context.basisFunctionType(), grid)

%}

