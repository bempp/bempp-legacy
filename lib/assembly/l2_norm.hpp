#ifndef bempp_l2_norm
#define bempp_l2_norm

#include "../common/common.hpp"
#include "../common/scalar_traits.hpp"

#include "evaluation_options.hpp"
#include "../fiber/quadrature_strategy.hpp"

namespace Fiber
{

template <typename ValueType> class Function;

} // namespace Fiber

namespace Bempp
{

class EvaluationOptions;
class GeometryFactory;
template <typename BasisFunctionType, typename ResultType> class GridFunction;

/** \brief Calculate the \f$L^2\f$-norm of the difference between a
 *  GridFunction and a Function.
 *
 *  This function can be used to estimate the absolute error of a numerical
 *  solution \p gridFunction against a known analytical solution \p refFunction.
 *
 *  The quadrature strategy \p quadStrategy is used to evaluate any
 *  necessary integrals; the \p options object controls the level of
 *  parallelism.
 */
template <typename BasisFunctionType, typename ResultType>
typename ScalarTraits<BasisFunctionType>::RealType L2NormOfDifference(
        const GridFunction<BasisFunctionType, ResultType>& gridFunction,
        const Fiber::Function<ResultType>& refFunction,
        const Fiber::QuadratureStrategy<
            BasisFunctionType, ResultType, GeometryFactory>& quadStrategy,
        const EvaluationOptions& options = EvaluationOptions());

/** \brief Calculate the absolute and relative \f$L^2\f$ errors of a solution.
 *
 *  This function calculates the absolute and relative \f$L^2\f$ norms
 *  of the difference between a GridFunction (typically a numerical
 *  solution of an integral equation) and a Function (typically a
 *  known analytical solution.
 *
 *  The quadrature strategy \p quadStrategy is used to evaluate any
 *  necessary integrals; the \p options object controls the level of
 *  parallelism. The results are returned in the output arguments \p
 *  absError and \p relError. */
template <typename BasisFunctionType, typename ResultType>
void estimateL2Error(
        const GridFunction<BasisFunctionType, ResultType>& gridFunction,
        const Fiber::Function<ResultType>& refFunction,
        const Fiber::QuadratureStrategy<
            BasisFunctionType, ResultType, GeometryFactory>& quadStrategy,
        const EvaluationOptions& options,
        typename ScalarTraits<BasisFunctionType>::RealType& absError,
        typename ScalarTraits<BasisFunctionType>::RealType& relError);

/** \brief Calculate the absolute and relative \f$L^2\f$ errors of a solution.
 *
 *  This is an overloaded version, provided for convenience. It is
 *  equivalent to the six-parameter version with \p options set to
 *  <tt>EvaluationOptions()</tt>. */
template <typename BasisFunctionType, typename ResultType>
void estimateL2Error(
        const GridFunction<BasisFunctionType, ResultType>& gridFunction,
        const Fiber::Function<ResultType>& refFunction,
        const Fiber::QuadratureStrategy<
            BasisFunctionType, ResultType, GeometryFactory>& quadStrategy,
        typename ScalarTraits<BasisFunctionType>::RealType& absError,
        typename ScalarTraits<BasisFunctionType>::RealType& relError);

} // namespace Bempp

#endif
