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

#ifndef fiber_quadrature_descriptor_selector_for_integral_operators_hpp
#define fiber_quadrature_descriptor_selector_for_integral_operators_hpp

#include "../common/common.hpp"

#include "double_quadrature_descriptor.hpp"

namespace Fiber {

/** \ingroup quadrature
 *  \brief Quadrature descriptor selector used during the
 *  discretization of boundary integral operators. */
template <typename CoordinateType>
class QuadratureDescriptorSelectorForIntegralOperators {
public:
  /** \brief Destructor */
  virtual ~QuadratureDescriptorSelectorForIntegralOperators() {}

  /** \brief Return the descriptor of the quadrature rule to be used
   *  for a particular pair of elements.
   *
   *  \param[in] testElementIndex Index of the test element.
   *  \param[in] trialElementIndex Index of the trial element.
   *  \param[in] nominalDistance
   *    This parameter is only relevant for quadrature selectors
   *    that vary the degree of accuracy of regular quadrature rules
   *    with interelement distance. If \p nominalDistance is
   *    non-negative, the quadrature rule should be chosen as if the
   *    distance between the centers of the test and trial element
   *    were equal to \p nominalDistance. Otherwise the selector
   *    should try to estimate the distance between the elements on
   *    its own and adjust the quadrature order accordingly. */
  virtual DoubleQuadratureDescriptor
  quadratureDescriptor(int testElementIndex, int trialElementIndex,
                       CoordinateType nominalDistance) const = 0;
};

} // namespace Fiber

#endif
