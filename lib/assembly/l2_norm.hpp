#ifndef bempp_l2_norm
#define bempp_l2_norm

#include "../common/common.hpp"
#include "../common/scalar_traits.hpp"

#include "../fiber/quadrature_strategy.hpp"

namespace Fiber
{

template <typename ValueType> class Function;
template <typename ValueType> class Function;

} // namespace Fiber

namespace Bempp
{

class EvaluationOptions;
class GeometryFactory;
template <typename BasisFunctionType, typename ResultType> class GridFunction;

template <typename BasisFunctionType, typename ResultType>
typename ScalarTraits<BasisFunctionType>::RealType L2NormOfDifference(
        const GridFunction<BasisFunctionType, ResultType>& gridFunction,
        const Fiber::Function<ResultType>& refFunction,
        const Fiber::QuadratureStrategy<
            BasisFunctionType, ResultType, GeometryFactory>& quadStrategy,
        const EvaluationOptions& options);

} // namespace Bempp

#endif
