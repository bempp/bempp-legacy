%{
#include "space/piecewise_constant_dual_mesh_scalar_space_barycentric.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs") piecewiseConstantDualMeshScalarSpaceBarycentric;
}

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
boost::shared_ptr<Space<BasisFunctionType> >
piecewiseConstantDualMeshScalarSpaceBarycentric(const boost::shared_ptr<const Grid>& grid,
        const GridSegment* segment = 0)
{
    typedef PiecewiseConstantDualMeshScalarSpaceBarycentric<BasisFunctionType> Type;
    if (segment)
        return boost::shared_ptr<Type>(new Type(grid, *segment));
    else
        return boost::shared_ptr<Type>(new Type(grid));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(piecewiseConstantDualMeshScalarSpaceBarycentric);
}

