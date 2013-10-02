%{
#include "space/piecewise_constant_dual_grid_scalar_space.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs") piecewiseConstantDualGridScalarSpace;
}

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
boost::shared_ptr<Space<BasisFunctionType> >
piecewiseConstantDualGridScalarSpace(const boost::shared_ptr<const Grid>& grid)
{
    typedef PiecewiseConstantDualGridScalarSpace<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(grid));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(piecewiseConstantDualGridScalarSpace);
}

