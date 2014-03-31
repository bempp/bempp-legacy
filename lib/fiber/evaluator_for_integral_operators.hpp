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

#ifndef fiber_evaluator_for_integral_operators_hpp
#define fiber_evaluator_for_integral_operators_hpp

#include "../common/armadillo_fwd.hpp"
#include "scalar_traits.hpp"

#include "../common/common.hpp"

namespace Fiber
{

/** \brief Interface of classes used to evaluate potentials defined by the
 *  application of integral operators to functions defined on surfaces. */
template <typename ResultType>
class EvaluatorForIntegralOperators
{
public:
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    enum Region {
        NEAR_FIELD, FAR_FIELD
    };

    virtual ~EvaluatorForIntegralOperators() {}

    /** \brief Evaluate the operator at a given list of points.
     *
     *  \param[in] region
     *   Obsolete. Always use FAR_FIELD.
     *
     *  \param[in] points
     *   2D array whose (i, j)th element is the ith coordinate of the jth point
     *   at which the potential should be evaluated.
     *
     *  \param[out] result
     *   On output, a 2D array whose (i, j)th element is the ith component of the
     *   potential at the jth point.
     */
    virtual void evaluate(Region region,
                          const arma::Mat<CoordinateType>& points,
                          arma::Mat<ResultType>& result) const = 0;
};

} // namespace Fiber

#endif
