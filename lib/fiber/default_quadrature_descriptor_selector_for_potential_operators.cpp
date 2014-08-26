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

#include "default_quadrature_descriptor_selector_for_potential_operators.hpp"

#include "default_local_assembler_for_operators_on_surfaces_utilities.hpp"
#include "explicit_instantiation.hpp"
#include "quadrature_options.hpp"
#include "raw_grid_geometry.hpp"
#include "shapeset.hpp"

namespace Fiber {

template <typename BasisFunctionType>
DefaultQuadratureDescriptorSelectorForPotentialOperators<BasisFunctionType>::
    DefaultQuadratureDescriptorSelectorForPotentialOperators(
        const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType> *>> &trialShapesets,
        const AccuracyOptionsEx &accuracyOptions)
    : m_rawGeometry(rawGeometry), m_trialShapesets(trialShapesets),
      m_accuracyOptions(accuracyOptions) {
  Utilities::checkConsistencyOfGeometryAndShapesets(*rawGeometry,
                                                    *trialShapesets);
  precalculateElementSizesAndCenters();
}

template <typename BasisFunctionType>
void DefaultQuadratureDescriptorSelectorForPotentialOperators<
    BasisFunctionType>::precalculateElementSizesAndCenters() {
  Utilities::precalculateElementSizesAndCentersForSingleGrid(
      *m_rawGeometry, m_elementSizesSquared, m_elementCenters,
      m_averageElementSize);
}

template <typename BasisFunctionType>
SingleQuadratureDescriptor
DefaultQuadratureDescriptorSelectorForPotentialOperators<BasisFunctionType>::
    quadratureDescriptor(const arma::Col<CoordinateType> &point,
                         int trialElementIndex,
                         CoordinateType nominalDistance) const {
  SingleQuadratureDescriptor desc;
  desc.vertexCount = m_rawGeometry->elementCornerCount(trialElementIndex);
  desc.order = order(point, trialElementIndex, nominalDistance);
  return desc;
}

template <typename BasisFunctionType>
int DefaultQuadratureDescriptorSelectorForPotentialOperators<
    BasisFunctionType>::order(const arma::Col<CoordinateType> &point,
                              int trialElementIndex,
                              CoordinateType nominalDistance) const {
  // Order required for exact quadrature on affine elements with a constant
  // kernel
  int trialBasisOrder = (*m_trialShapesets)[trialElementIndex]->order();
  int defaultQuadOrder = 2 * trialBasisOrder;

  CoordinateType normalisedDistance;
  if (nominalDistance < 0.) {
    CoordinateType elementSizeSquared =
        m_elementSizesSquared[trialElementIndex];
    CoordinateType distanceSquared =
        pointElementDistanceSquared(point, trialElementIndex);
    CoordinateType normalisedDistanceSquared =
        distanceSquared / elementSizeSquared;
    normalisedDistance = sqrt(normalisedDistanceSquared);
  } else
    normalisedDistance = nominalDistance / m_averageElementSize;

  const QuadratureOptions &options =
      m_accuracyOptions.singleRegular(normalisedDistance);
  return options.quadratureOrder(defaultQuadOrder);
}

template <typename BasisFunctionType>
inline typename DefaultQuadratureDescriptorSelectorForPotentialOperators<
    BasisFunctionType>::CoordinateType
DefaultQuadratureDescriptorSelectorForPotentialOperators<BasisFunctionType>::
    pointElementDistanceSquared(const arma::Col<CoordinateType> &point,
                                int trialElementIndex) const {
  // TODO: optimize
  arma::Col<CoordinateType> diff =
      point - m_elementCenters.col(trialElementIndex);
  return arma::dot(diff, diff);
}

template <typename BasisFunctionType>
SingleQuadratureDescriptor
DefaultQuadratureDescriptorSelectorForPotentialOperators<BasisFunctionType>::
    farFieldQuadratureDescriptor(
        const Shapeset<BasisFunctionType> &trialShapeset,
        int trialElementCornerCount) const {
  SingleQuadratureDescriptor desc;
  desc.vertexCount = trialElementCornerCount;
  int elementOrder = trialShapeset.order();
  // Order required for exact quadrature on affine elements with kernel
  // approximated by a polynomial of order identical with that of the shapeset
  int defaultQuadratureOrder = 2 * elementOrder;
  desc.order =
      m_accuracyOptions.singleRegular().quadratureOrder(defaultQuadratureOrder);
  return desc;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(
    DefaultQuadratureDescriptorSelectorForPotentialOperators);

} // namespace Fiber
