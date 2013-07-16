%{
#include "space/piecewise_linear_discontinuous_scalar_space.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs") piecewiseLinearDiscontinuousScalarSpace;
}

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
boost::shared_ptr<Space<BasisFunctionType> >
piecewiseLinearDiscontinuousScalarSpace(const boost::shared_ptr<const Grid>& grid,
    const GridSegment* segment = 0)
{
    typedef PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType> Type;
    if (segment)
        return boost::shared_ptr<Type>(new Type(grid, *segment));
    else
        return boost::shared_ptr<Type>(new Type(grid));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(piecewiseLinearDiscontinuousScalarSpace);
}

