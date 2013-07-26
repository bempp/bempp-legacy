%{
#include "space/raviart_thomas_0_vector_space.hpp"
%}

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
boost::shared_ptr<Space<BasisFunctionType> >
raviartThomas0VectorSpace(const boost::shared_ptr<const Grid>& grid,
                          const GridSegment* segment = 0,
                          bool putDofsOnBoundaries = false,
                          int dofMode = Bempp::EDGE_ON_SEGMENT)
{
    typedef RaviartThomas0VectorSpace<BasisFunctionType> Type;
    if (segment)
        return boost::shared_ptr<Type>(new Type(grid, *segment, 
                                                putDofsOnBoundaries,
                                                dofMode));
    else
        return boost::shared_ptr<Type>(new Type(grid, putDofsOnBoundaries));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(raviartThomas0VectorSpace);
}

