%{
#include "space/piecewise_polynomial_discontinuous_scalar_space.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs") piecewisePolynomialDiscontinuousScalarSpace;
}

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
boost::shared_ptr<Space<BasisFunctionType> >
piecewisePolynomialDiscontinuousScalarSpace(
    const boost::shared_ptr<const Grid>& grid,
    int polynomialOrder,
    const GridSegment* segment = 0,
    bool strictlyOnSegment = false)
{
    typedef PiecewisePolynomialDiscontinuousScalarSpace<BasisFunctionType> Type;
    if (segment)
        return boost::shared_ptr<Type>(new Type(grid, polynomialOrder, *segment,
                                                strictlyOnSegment));
    else
        return boost::shared_ptr<Type>(new Type(grid, polynomialOrder));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(piecewisePolynomialDiscontinuousScalarSpace);
}

