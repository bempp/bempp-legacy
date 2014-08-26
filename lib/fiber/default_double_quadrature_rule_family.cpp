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

#include "default_double_quadrature_rule_family.hpp"

#include "explicit_instantiation.hpp"
#include "numerical_quadrature.hpp"

namespace Fiber {

template <typename CoordinateType>
void DefaultDoubleQuadratureRuleFamily<CoordinateType>::
    fillQuadraturePointsAndWeights(const DoubleQuadratureDescriptor &desc,
                                   arma::Mat<CoordinateType> &testPoints,
                                   arma::Mat<CoordinateType> &trialPoints,
                                   std::vector<CoordinateType> &testWeights,
                                   std::vector<CoordinateType> &trialWeights,
                                   bool &isTensor) const {
  const ElementPairTopology &topology = desc.topology;
  if (topology.type == ElementPairTopology::Disjoint) {
    // Create a tensor rule
    fillSingleQuadraturePointsAndWeights(
        topology.testVertexCount, desc.testOrder, testPoints, testWeights);
    fillSingleQuadraturePointsAndWeights(
        topology.trialVertexCount, desc.trialOrder, trialPoints, trialWeights);
    isTensor = true;
  } else {
    // Create a non-tensor rule, leaving trialWeights empty
    fillDoubleSingularQuadraturePointsAndWeights(desc, testPoints, trialPoints,
                                                 testWeights);
    trialWeights.clear();
    isTensor = false;
  }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT_REAL_ONLY(
    DefaultDoubleQuadratureRuleFamily);

} // namespace Fiber
