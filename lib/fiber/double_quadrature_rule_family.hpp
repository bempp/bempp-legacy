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

#ifndef fiber_double_quadrature_rule_family_hpp
#define fiber_double_quadrature_rule_family_hpp

#include "../common/common.hpp"
#include "../common/armadillo_fwd.hpp"
#include "double_quadrature_descriptor.hpp"

#include <vector>

namespace Fiber {

/** \brief Family of quadrature rules over pairs of elements.
 *
 *  This class represents a family of quadrature rules that allows
 *  integration of possibly singular functions f(x, y) over pairs of
 *  reference elements (triangles or squares).  Individual quadrature
 *  rules may be either separable (tensor-product), i.e. of the form
 *  \[f
 *  \int_S \int_T f(x, y) \, \mathrm d S(x) \, \mathrm d T(x) \approx
 *  \sum_{i=1}^M \sum_{j=1}^N f(x_i, y_j) w_i w_j,
 *  \f]
 *  or nonseparable (non-tensor-product), i.e. of the form
 *  \[f
 *  \int_S \int_T f(x, y) \, \mathrm d S(x) \, \mathrm d T(x) \approx
 *  \sum_{k=1}^K f(x_k, y_k) w_k.
 *  \f]
 *  Above, each of \f$S\f$ and \f$T\f$ denotes either the unit
 *  triangle or the unit square.
 */
template <typename CoordinateType> class DoubleQuadratureRuleFamily {
public:
  /** \brief Destructor. */
  virtual ~DoubleQuadratureRuleFamily() {}

  /** \brief Fill arrays of quadrature points and weights.
   *
   *  \param[in]  desc    Quadrature descriptor.
   *  \param[out] testPoints
   *    2D array whose (i, j)th element will contain the ith (local)
   *    coordinate of the jth quadrature point on the test element.
   *  \param[out] trialPoints
   *    2D array whose (i, j)th element will contain the ith (local)
   *    coordinate of the jth quadrature point on the trial element.
   *  \param[out] testWeights
   *    If \p isTensor is set to \c true, this parameter will be a
   *    vector of test point weights of a separable quadrature rule;
   *    otherwise it will be a vector of point weights of a
   *    nonseparable quadarture rule.
   *  \param[out] trialWeights
   *    If \p isTensor is set to \c true, this parameter will be a
   *    vector of trial point weights of a separable quadrature rule;
   *    otherwise it will be empty.
   *  \param[out] isTensor
   *    If set to \c true, the returned quadrature rule is a
   *    tensor-product rule, otherwise it is not. */
  virtual void
  fillQuadraturePointsAndWeights(const DoubleQuadratureDescriptor &desc,
                                 arma::Mat<CoordinateType> &testPoints,
                                 arma::Mat<CoordinateType> &trialPoints,
                                 std::vector<CoordinateType> &testWeights,
                                 std::vector<CoordinateType> &trialWeights,
                                 bool &isTensor) const = 0;
};

} // namespace Fiber

#endif
