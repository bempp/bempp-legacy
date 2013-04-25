%{
#include "space/piecewise_linear_discontinuous_scalar_space.hpp"
%}

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
boost::shared_ptr<Space<BasisFunctionType> >
piecewiseLinearDiscontinuousScalarSpace(const boost::shared_ptr<const Grid>& grid)
{
    typedef PiecewiseLinearDiscontinuousScalarSpace<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(grid));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(piecewiseLinearDiscontinuousScalarSpace);
}

