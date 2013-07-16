%{
#include "space/piecewise_polynomial_continuous_scalar_space.hpp"
%}

namespace Bempp
{
%feature("compactdefaultargs") piecewisePolynomialContinuousScalarSpace;
}

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
boost::shared_ptr<Space<BasisFunctionType> >
piecewisePolynomialContinuousScalarSpace(
    const boost::shared_ptr<const Grid>& grid,
    int polynomialOrder,
    const GridSegment* segment = 0)
{
    typedef PiecewisePolynomialContinuousScalarSpace<BasisFunctionType> Type;
    if (segment)
        return boost::shared_ptr<Type>(new Type(grid, polynomialOrder, *segment));
    else
        return boost::shared_ptr<Type>(new Type(grid, polynomialOrder));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(piecewisePolynomialContinuousScalarSpace);
}

