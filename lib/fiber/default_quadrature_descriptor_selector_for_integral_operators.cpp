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

#include "default_quadrature_descriptor_selector_for_integral_operators.hpp"

#include "default_local_assembler_for_operators_on_surfaces_utilities.hpp"
#include "explicit_instantiation.hpp"
#include "quadrature_options.hpp"
#include "raw_grid_geometry.hpp"
#include "shapeset.hpp"

namespace Fiber
{

template <typename BasisFunctionType>
DefaultQuadratureDescriptorSelectorForIntegralOperators<BasisFunctionType>::
DefaultQuadratureDescriptorSelectorForIntegralOperators(
        const shared_ptr<const RawGridGeometry<CoordinateType> >& testRawGeometry,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& trialRawGeometry,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType>*> >& testShapesets,
        const shared_ptr<const std::vector<
            const Shapeset<BasisFunctionType>*> >& trialShapesets,
        const AccuracyOptionsEx& accuracyOptions) :
    m_testRawGeometry(testRawGeometry),
    m_trialRawGeometry(trialRawGeometry),
    m_testShapesets(testShapesets),
    m_trialShapesets(trialShapesets),
    m_accuracyOptions(accuracyOptions)
{
    Utilities::checkConsistencyOfGeometryAndShapesets(
        *testRawGeometry, *testShapesets);
    Utilities::checkConsistencyOfGeometryAndShapesets(
        *trialRawGeometry, *trialShapesets);
    precalculateElementSizesAndCenters();
}

template <typename BasisFunctionType>
inline bool
DefaultQuadratureDescriptorSelectorForIntegralOperators<BasisFunctionType>::
testAndTrialGridsAreIdentical() const
{
    return m_testRawGeometry.get() == m_trialRawGeometry.get();
}

template <typename BasisFunctionType>
void
DefaultQuadratureDescriptorSelectorForIntegralOperators<BasisFunctionType>::
precalculateElementSizesAndCenters()
{
    CoordinateType averageTestElementSize;
    Utilities::precalculateElementSizesAndCentersForSingleGrid(
                *m_testRawGeometry,
                m_testElementSizesSquared, m_testElementCenters,
                averageTestElementSize);
    if (testAndTrialGridsAreIdentical()) {
        m_trialElementSizesSquared = m_testElementSizesSquared;
        m_trialElementCenters = m_testElementCenters;
        m_averageElementSize = averageTestElementSize;
    } else {
        CoordinateType averageTrialElementSize;
        Utilities::precalculateElementSizesAndCentersForSingleGrid(
                    *m_trialRawGeometry,
                    m_trialElementSizesSquared, m_trialElementCenters,
                    averageTrialElementSize);
        m_averageElementSize =
                (averageTestElementSize + averageTrialElementSize) / 2.;
    }
}

template <typename BasisFunctionType>
DoubleQuadratureDescriptor
DefaultQuadratureDescriptorSelectorForIntegralOperators<BasisFunctionType>::
quadratureDescriptor(
        int testElementIndex, int trialElementIndex,
        CoordinateType nominalDistance) const
{
    DoubleQuadratureDescriptor desc;

    // Get corner indices of the specified elements
    arma::Col<int> testElementCornerIndices =
            m_testRawGeometry->elementCornerIndices(testElementIndex);
    arma::Col<int> trialElementCornerIndices =
            m_trialRawGeometry->elementCornerIndices(trialElementIndex);
    if (testAndTrialGridsAreIdentical()) {
        desc.topology = determineElementPairTopologyIn3D(
                    testElementCornerIndices, trialElementCornerIndices);
    } else {
        desc.topology.testVertexCount = testElementCornerIndices.n_rows;
        desc.topology.trialVertexCount = trialElementCornerIndices.n_rows;
        desc.topology.type = ElementPairTopology::Disjoint;
    }

    if (desc.topology.type == ElementPairTopology::Disjoint) {
        getRegularOrders(testElementIndex, trialElementIndex,
                         desc.testOrder, desc.trialOrder,
                         nominalDistance);
    } else { // singular integral
        desc.testOrder = singularOrder(testElementIndex, TEST);
        desc.trialOrder = singularOrder(trialElementIndex, TRIAL);
    }

    return desc;
}

template <typename BasisFunctionType>
void
DefaultQuadratureDescriptorSelectorForIntegralOperators<BasisFunctionType>::
getRegularOrders(int testElementIndex, int trialElementIndex,
                 int& testQuadOrder, int& trialQuadOrder,
                 CoordinateType nominalDistance) const
{
    // TODO:
    // 1. Check the size of elements and the distance between them
    //    and estimate the variability of the kernel
    // 2. Take into account the fact that elements might be isoparametric.

    // Order required for exact quadrature on affine elements with a constant kernel
    int testBasisOrder = (*m_testShapesets)[testElementIndex]->order();
    int trialBasisOrder = (*m_trialShapesets)[trialElementIndex]->order();
    testQuadOrder = testBasisOrder;
    trialQuadOrder = trialBasisOrder;

    CoordinateType normalisedDistance;
    if (nominalDistance < 0.) {
        CoordinateType testElementSizeSquared =
                m_testElementSizesSquared[testElementIndex];
        CoordinateType trialElementSizeSquared =
                m_trialElementSizesSquared[trialElementIndex];
        CoordinateType distanceSquared =
                elementDistanceSquared(testElementIndex, trialElementIndex);
        CoordinateType normalisedDistanceSquared =
                distanceSquared / std::max(testElementSizeSquared,
                                           trialElementSizeSquared);
        normalisedDistance = sqrt(normalisedDistanceSquared);
    } else
        normalisedDistance = nominalDistance / m_averageElementSize;

    const QuadratureOptions& options =
            m_accuracyOptions.doubleRegular(normalisedDistance);
    testQuadOrder = options.quadratureOrder(testQuadOrder);
    trialQuadOrder = options.quadratureOrder(trialQuadOrder);
}

template <typename BasisFunctionType>
int
DefaultQuadratureDescriptorSelectorForIntegralOperators<BasisFunctionType>::
singularOrder(int elementIndex, ElementType elementType) const
{
    // TODO:
    // 1. Check the size of elements and estimate the variability of the
    //    (Sauter-Schwab-transformed) kernel
    // 2. Take into account the fact that elements might be isoparametric.

    const QuadratureOptions& options = m_accuracyOptions.doubleSingular();

    int elementOrder = (elementType == TEST ?
                            (*m_testShapesets)[elementIndex]->order() :
                            (*m_trialShapesets)[elementIndex]->order());
    int defaultAccuracyOrder = elementOrder + 5;
    return options.quadratureOrder(defaultAccuracyOrder);
}

template <typename BasisFunctionType>
inline
typename DefaultQuadratureDescriptorSelectorForIntegralOperators<
BasisFunctionType>::CoordinateType
DefaultQuadratureDescriptorSelectorForIntegralOperators<BasisFunctionType>::
elementDistanceSquared(
        int testElementIndex, int trialElementIndex) const
{
    CoordinateType result = 0.;
    const int dimWorld = 3;
    for (int d = 0; d < dimWorld; ++d) {
        CoordinateType diff = m_trialElementCenters(d, trialElementIndex) -
                m_testElementCenters(d, testElementIndex);
        result += diff * diff;
    }
    return result;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(
    DefaultQuadratureDescriptorSelectorForIntegralOperators);

} // namespace Fiber
