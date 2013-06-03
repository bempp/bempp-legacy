%{
#include "space/piecewise_polynomial_continuous_scalar_space.hpp"
%}

%inline %{
namespace Bempp
{

template <typename BasisFunctionType>
boost::shared_ptr<Space<BasisFunctionType> >
piecewisePolynomialContinuousScalarSpace(
    const boost::shared_ptr<const Grid>& grid,
    int polynomialOrder)
{
    typedef PiecewisePolynomialContinuousScalarSpace<BasisFunctionType> Type;
    return boost::shared_ptr<Type>(new Type(grid, polynomialOrder));
}

} // namespace Bempp
%}

namespace Bempp
{
BEMPP_INSTANTIATE_SYMBOL_TEMPLATED_ON_BASIS(piecewisePolynomialContinuousScalarSpace);
}

