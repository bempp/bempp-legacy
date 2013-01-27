%{
#include "space/raviart_thomas_0_vector_space.hpp"
%}

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
boost::shared_ptr<Space<BasisFunctionType> >
raviartThomas0VectorSpace(const boost::shared_ptr<const Grid>& grid)
{
    typedef RaviartThomas0VectorSpace<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(grid));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(raviartThomas0VectorSpace);
}

