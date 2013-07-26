%{
#include "space/piecewise_linear_continuous_scalar_space.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs") piecewiseLinearContinuousScalarSpace;
}

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
boost::shared_ptr<Space<BasisFunctionType> >
piecewiseLinearContinuousScalarSpace(
    const boost::shared_ptr<const Grid>& grid,
    const GridSegment* segment = 0,
    bool strictlyOnSegment = false)
{
    typedef PiecewiseLinearContinuousScalarSpace<BasisFunctionType> Type;
    if (segment)
        return boost::shared_ptr<Type>(new Type(grid, *segment, 
                                                strictlyOnSegment));
    else
        return boost::shared_ptr<Type>(new Type(grid));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(piecewiseLinearContinuousScalarSpace);
}

