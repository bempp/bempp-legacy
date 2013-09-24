%{
#include "space/piecewise_constant_dual_grid_discontinuous_scalar_space.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs") piecewiseConstantDualGridDiscontinuousScalarSpace;
}

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
boost::shared_ptr<Space<BasisFunctionType> >
piecewiseConstantDualGridDiscontinuousScalarSpace(const boost::shared_ptr<const Grid>& grid)
{
    typedef PiecewiseConstantDualGridDiscontinuousScalarSpace<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(grid));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(piecewiseConstantDualGridDiscontinuousScalarSpace);
}

