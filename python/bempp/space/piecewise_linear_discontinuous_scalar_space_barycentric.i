%{
#include "space/piecewise_linear_discontinuous_scalar_space_barycentric.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs") piecewiseLinearDisontinuousScalarSpaceBarycentric;
}

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
boost::shared_ptr<Space<BasisFunctionType> >
piecewiseLinearDiscontinuousScalarSpaceBarycentric(
    const boost::shared_ptr<const Grid>& grid,
    const GridSegment* segment = 0,
    bool strictlyOnSegment = false)
{
    typedef PiecewiseLinearDiscontinuousScalarSpaceBarycentric<BasisFunctionType> Type;
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
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(piecewiseLinearDiscontinuousScalarSpaceBarycentric);
}
