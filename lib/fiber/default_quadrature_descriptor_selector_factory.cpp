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

#include "default_quadrature_descriptor_selector_factory.hpp"

#include "default_quadrature_descriptor_selector_for_grid_functions.hpp"
#include "default_quadrature_descriptor_selector_for_integral_operators.hpp"
#include "default_quadrature_descriptor_selector_for_local_operators.hpp"
#include "default_quadrature_descriptor_selector_for_potential_operators.hpp"
#include "explicit_instantiation.hpp"

#include "../common/boost_make_shared_fwd.hpp"

namespace Fiber {

template <typename BasisFunctionType>
DefaultQuadratureDescriptorSelectorFactory<BasisFunctionType>::
    DefaultQuadratureDescriptorSelectorFactory(
        const AccuracyOptionsEx &accuracyOptions)
    : m_accuracyOptions(accuracyOptions) {}

template <typename BasisFunctionType>
shared_ptr<QuadratureDescriptorSelectorForGridFunctions<
    typename DefaultQuadratureDescriptorSelectorFactory<
        BasisFunctionType>::CoordinateType>>
DefaultQuadratureDescriptorSelectorFactory<BasisFunctionType>::
    makeQuadratureDescriptorSelectorForGridFunctions(
        const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
            &testShapesets) const {
  typedef DefaultQuadratureDescriptorSelectorForGridFunctions<BasisFunctionType>
      Selector;
  return boost::make_shared<Selector>(rawGeometry, testShapesets,
                                      m_accuracyOptions);
}

template <typename BasisFunctionType>
shared_ptr<QuadratureDescriptorSelectorForIntegralOperators<
    typename DefaultQuadratureDescriptorSelectorFactory<
        BasisFunctionType>::CoordinateType>>
DefaultQuadratureDescriptorSelectorFactory<BasisFunctionType>::
    makeQuadratureDescriptorSelectorForIntegralOperators(
        const shared_ptr<const RawGridGeometry<CoordinateType>>
            &testRawGeometry,
        const shared_ptr<const RawGridGeometry<CoordinateType>>
            &trialRawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
            &testShapesets,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
            &trialShapesets) const {
  typedef DefaultQuadratureDescriptorSelectorForIntegralOperators<
      BasisFunctionType> Selector;
  return boost::make_shared<Selector>(testRawGeometry, trialRawGeometry,
                                      testShapesets, trialShapesets,
                                      m_accuracyOptions);
}

template <typename BasisFunctionType>
shared_ptr<QuadratureDescriptorSelectorForLocalOperators<
    typename DefaultQuadratureDescriptorSelectorFactory<
        BasisFunctionType>::CoordinateType>>
DefaultQuadratureDescriptorSelectorFactory<BasisFunctionType>::
    makeQuadratureDescriptorSelectorForLocalOperators(
        const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
            &testShapesets,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
            &trialShapesets) const {
  typedef DefaultQuadratureDescriptorSelectorForLocalOperators<
      BasisFunctionType> Selector;
  return boost::make_shared<Selector>(rawGeometry, testShapesets,
                                      trialShapesets, m_accuracyOptions);
}

template <typename BasisFunctionType>
shared_ptr<QuadratureDescriptorSelectorForPotentialOperators<BasisFunctionType>>
DefaultQuadratureDescriptorSelectorFactory<BasisFunctionType>::
    makeQuadratureDescriptorSelectorForPotentialOperators(
        const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
            &trialShapesets) const {
  typedef DefaultQuadratureDescriptorSelectorForPotentialOperators<
      BasisFunctionType> Selector;
  return boost::make_shared<Selector>(rawGeometry, trialShapesets,
                                      m_accuracyOptions);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(
    DefaultQuadratureDescriptorSelectorFactory);

} // namespace Fiber
