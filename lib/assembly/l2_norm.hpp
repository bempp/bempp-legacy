// Copyright (C) 2011-2012 by the BEM++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef bempp_l2_norm
#define bempp_l2_norm

/** \file l2_norm.hpp Calculation of \f$L^2\f$-norms of grid functions. */

#include "../common/common.hpp"
#include "../common/scalar_traits.hpp"

#include "evaluation_options.hpp"
#include "../fiber/quadrature_strategy.hpp"

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename ValueType> class Function;
/** \endcond */

} // namespace Fiber

namespace Bempp {

/** \cond FORWARD_DECL */
class EvaluationOptions;
class GeometryFactory;
template <typename BasisFunctionType, typename ResultType> class GridFunction;
/** \endcond */

/** \relates GridFunction
 *  \brief Calculate the \f$L^2\f$-norm of the difference between a
 *  GridFunction and a Function.
 *
 *  This function can be used to estimate the absolute error of a numerical
 *  solution \p gridFunction against a known analytical solution \p refFunction.
 *
 *  The quadrature strategy \p quadStrategy is used to evaluate any
 *  necessary integrals; the \p options object controls the level of
 *  parallelism.
 *
 *  \note If you use the numerical quadrature strategy, you may need to increase
 *  the quadrature order for regular integrals on single elements by at least 2
 *  to ensure that this function produces results with sufficient accuracy.
 */
template <typename BasisFunctionType, typename ResultType>
typename ScalarTraits<BasisFunctionType>::RealType L2NormOfDifference(
    const GridFunction<BasisFunctionType, ResultType> &gridFunction,
    const Fiber::Function<ResultType> &refFunction,
    const Fiber::QuadratureStrategy<BasisFunctionType, ResultType,
                                    GeometryFactory> &quadStrategy,
    const EvaluationOptions &options = EvaluationOptions());

/** \relates GridFunction
 *  \brief Calculate the absolute and relative \f$L^2\f$ errors of a solution.
 *
 *  This function calculates the absolute and relative \f$L^2\f$ norms
 *  of the difference between a GridFunction (typically a numerical
 *  solution of an integral equation) and a Function (typically a
 *  known analytical solution.
 *
 *  The quadrature strategy \p quadStrategy is used to evaluate any
 *  necessary integrals; the \p options object controls the level of
 *  parallelism. The results are returned in the output arguments \p
 *  absError and \p relError.
 *
 *  \note If you use the numerical quadrature strategy, you may need to increase
 *  the quadrature order for regular integrals on single elements by at least 2
 *  to ensure that this function produces results with sufficient accuracy.
 */
template <typename BasisFunctionType, typename ResultType>
void
estimateL2Error(const GridFunction<BasisFunctionType, ResultType> &gridFunction,
                const Fiber::Function<ResultType> &refFunction,
                const Fiber::QuadratureStrategy<BasisFunctionType, ResultType,
                                                GeometryFactory> &quadStrategy,
                const EvaluationOptions &options,
                typename ScalarTraits<BasisFunctionType>::RealType &absError,
                typename ScalarTraits<BasisFunctionType>::RealType &relError);

/** \relates GridFunction
 *  \brief Calculate the absolute and relative \f$L^2\f$ errors of a solution.
 *
 *  This is an overloaded version, provided for convenience. It is
 *  equivalent to the six-parameter version with \p options set to
 *  <tt>EvaluationOptions()</tt>. */
template <typename BasisFunctionType, typename ResultType>
void
estimateL2Error(const GridFunction<BasisFunctionType, ResultType> &gridFunction,
                const Fiber::Function<ResultType> &refFunction,
                const Fiber::QuadratureStrategy<BasisFunctionType, ResultType,
                                                GeometryFactory> &quadStrategy,
                typename ScalarTraits<BasisFunctionType>::RealType &absError,
                typename ScalarTraits<BasisFunctionType>::RealType &relError);

} // namespace Bempp

#endif
