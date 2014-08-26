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

#include "default_quadrature_descriptor_selector_for_local_operators.hpp"

#include "default_local_assembler_for_operators_on_surfaces_utilities.hpp"
#include "explicit_instantiation.hpp"
#include "quadrature_options.hpp"
#include "raw_grid_geometry.hpp"
#include "shapeset.hpp"

namespace Fiber {

template <typename BasisFunctionType>
DefaultQuadratureDescriptorSelectorForLocalOperators<BasisFunctionType>::
    DefaultQuadratureDescriptorSelectorForLocalOperators(
        const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType> *>> &testShapesets,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType> *>> &trialShapesets,
        const AccuracyOptionsEx &accuracyOptions)
    : m_rawGeometry(rawGeometry), m_testShapesets(testShapesets),
      m_trialShapesets(trialShapesets), m_accuracyOptions(accuracyOptions) {
  Utilities::checkConsistencyOfGeometryAndShapesets(*rawGeometry,
                                                    *testShapesets);
  Utilities::checkConsistencyOfGeometryAndShapesets(*rawGeometry,
                                                    *trialShapesets);
}

template <typename BasisFunctionType>
SingleQuadratureDescriptor DefaultQuadratureDescriptorSelectorForLocalOperators<
    BasisFunctionType>::quadratureDescriptor(int elementIndex) const {
  SingleQuadratureDescriptor desc;

  // Get number of corners of the specified element
  desc.vertexCount = m_rawGeometry->elementCornerCount(elementIndex);

  // Determine integrand's order and required quadrature order
  const int expressionOrder = (*m_testShapesets)[elementIndex]->order() +
                              (*m_trialShapesets)[elementIndex]->order();
  desc.order = expressionOrder;

  return desc;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(
    DefaultQuadratureDescriptorSelectorForLocalOperators);

} // namespace Fiber
