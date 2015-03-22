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

#ifndef fiber_default_double_quadrature_rule_family_hpp
#define fiber_default_double_quadrature_rule_family_hpp

#include "double_quadrature_rule_family.hpp"
#include "types.hpp"

namespace Fiber {

/** \brief Default family of quadrature rules over pairs of elements.
 *
 *  Regular integrals are treated with tensor-product quadrature
 *  rules, singular integrals with non-tensor-product rules taken from
 *  Sauter and Schwab, "Boundary Element Methods", Springer 2011. */
template <typename CoordinateType>
class DefaultDoubleQuadratureRuleFamily
    : public DoubleQuadratureRuleFamily<CoordinateType> {
public:
  virtual ~DefaultDoubleQuadratureRuleFamily() {}

  virtual void
  fillQuadraturePointsAndWeights(const DoubleQuadratureDescriptor &desc,
                                 Matrix<CoordinateType> &testPoints,
                                 Matrix<CoordinateType> &trialPoints,
                                 std::vector<CoordinateType> &testWeights,
                                 std::vector<CoordinateType> &trialWeights,
                                 bool &isTensor) const;
};

} // namespace Fiber

#endif
