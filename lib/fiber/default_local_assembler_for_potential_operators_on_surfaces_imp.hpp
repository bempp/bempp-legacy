// Copyright (C) 2011-2012 by the Bem++ Authors
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

// Keep IDEs happy
#include "default_local_assembler_for_potential_operators_on_surfaces.hpp"

#include "basis.hpp"
#include "kernel_trial_integral.hpp"
#include "numerical_kernel_trial_integrator.hpp"
#include "raw_grid_geometry.hpp"
#include "serial_blas_region.hpp"

#include <cassert>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

#include "../common/auto_timer.hpp"

namespace Fiber
{

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces(
        const arma::Mat<CoordinateType>& points,
        const shared_ptr<const GeometryFactory>& geometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType> >& rawGeometry,
        const shared_ptr<const std::vector<const Basis<BasisFunctionType>*> >& trialBases,
        const shared_ptr<const CollectionOfKernels<KernelType> >& kernels,
        const shared_ptr<const CollectionOfBasisTransformations<CoordinateType> >& trialTransformations,
        const shared_ptr<const KernelTrialIntegral<BasisFunctionType, KernelType, ResultType> >& integral,
        const ParallelizationOptions& parallelizationOptions,
        VerbosityLevel::Level verbosityLevel,
        const AccuracyOptionsEx& accuracyOptions) :
    m_points(points),
    m_geometryFactory(geometryFactory),
    m_rawGeometry(rawGeometry),
    m_trialBases(trialBases),
    m_kernels(kernels),
    m_trialTransformations(trialTransformations),
    m_integral(integral),
    m_parallelizationOptions(parallelizationOptions),
    m_verbosityLevel(verbosityLevel),
    m_accuracyOptions(accuracyOptions)
{
    checkConsistencyOfGeometryAndBases(*rawGeometry, *trialBases);
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
~DefaultLocalAssemblerForPotentialOperatorsOnSurfaces()
{
    // Note: obviously the destructor is assumed to be called only after
    // all threads have ceased using the assembler!

    for (typename IntegratorMap::const_iterator it = m_kernelTrialIntegrators.begin();
         it != m_kernelTrialIntegrators.end(); ++it)
        delete it->second;
    m_kernelTrialIntegrators.clear();
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
checkConsistencyOfGeometryAndBases(
        const RawGridGeometry<CoordinateType>& rawGeometry,
        const std::vector<const Basis<BasisFunctionType>*>& bases) const
{
    if (rawGeometry.vertices().n_rows != 3)
        throw std::invalid_argument(
            "DefaultLocalAssemblerForPotentialOperatorsOnSurfaces::"
            "checkConsistencyOfGeometryAndBases(): "
            "vertex coordinates must be three-dimensional");
    const size_t elementCount = rawGeometry.elementCornerIndices().n_cols;
    if (rawGeometry.elementCornerIndices().n_rows < 3 ||
            4 < rawGeometry.elementCornerIndices().n_rows)
        throw std::invalid_argument(
            "DefaultLocalAssemblerForPotentialOperatorsOnSurfaces::"
            "checkConsistencyOfGeometryAndBases(): "
            "Elements must have either 3 or 4 corners");
    if (!rawGeometry.auxData().is_empty() &&
            rawGeometry.auxData().n_cols != elementCount)
        throw std::invalid_argument(
            "DefaultLocalAssemblerForPotentialOperatorsOnSurfaces::"
            "checkConsistencyOfGeometryAndBases(): "
            "number of columns of auxData must match that of "
            "elementCornerIndices");
    if (bases.size() != elementCount)
        throw std::invalid_argument(
            "DefaultLocalAssemblerForPotentialOperatorsOnSurfaces::"
            "checkConsistencyOfGeometryAndBases(): "
            "size of bases must match the number of columns of "
            "elementCornerIndices");
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
evaluateLocalContributions(
        const std::vector<int>& pointIndices,
        int trialElementIndex,
        LocalDofIndex localTrialDofIndex,
        std::vector<arma::Mat<ResultType> >& result,
        CoordinateType nominalDistance)
{
    typedef Basis<BasisFunctionType> Basis;

    const int pointCount = pointIndices.size();
    // const int componentCount = m_integral->resultDimension();

    // Get basis
    const Basis& trialBasis = *((*m_trialBases)[trialElementIndex]);
    assert(localTrialDofIndex == ALL_DOFS ||
           (localTrialDofIndex >= 0 && localTrialDofIndex < trialBasis.size()));
//    const Basis& trialBasis = *m_trialBases[trialElementIndex];
//    const int dofCount = (localTrialDofIndex == ALL_DOFS ?
//                              trialBasis.size() : 1);
//    assert(localTrialDofIndex == ALL_DOFS ||
//           (localTrialDofIndex >= 0 && localTrialDofIndex < dofCount));
    result.resize(pointCount);

    // Select integrators
    typedef const Integrator* QuadVariant;
    std::vector<QuadVariant> quadVariants(pointCount);
    for (int i = 0; i < pointCount; ++i) {
        const Integrator* integrator =
                &selectIntegrator(pointIndices[i], trialElementIndex,
                                  nominalDistance);
        quadVariants[i] = QuadVariant(integrator);
    }

    // Integration will proceed in batches of test elements having the same
    // "quadrature variant", i.e. integrator and basis

    // Find all the unique quadrature variants present
    typedef std::set<QuadVariant> QuadVariantSet;
    // Set of unique quadrature variants
    QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

    std::vector<int> activePointIndices;
    activePointIndices.reserve(pointCount);
    std::vector<arma::Mat<ResultType>*> activeLocalResults;
    activeLocalResults.reserve(pointCount);

    // Now loop over unique quadrature variants
    for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
         it != uniqueQuadVariants.end(); ++it) {
        const QuadVariant activeQuadVariant = *it;
        const Integrator& activeIntegrator = *activeQuadVariant;

        // Find all the test elements for which quadrature should proceed
        // according to the current quadrature variant
        activePointIndices.clear();
        activeLocalResults.clear();
        for (int point = 0; point < pointCount; ++point)
            if (quadVariants[point] == activeQuadVariant) {
                activePointIndices.push_back(pointIndices[point]);
                activeLocalResults.push_back(&result[point]);
            }

        // Integrate!
        activeIntegrator.integrate(activePointIndices, trialElementIndex,
                                   trialBasis, localTrialDofIndex,
                                   activeLocalResults);
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
evaluateLocalContributions(
        int pointIndex,
        int componentIndex,
        const std::vector<int>& trialElementIndices,
        std::vector<arma::Mat<ResultType> >& result,
        CoordinateType nominalDistance)
{
    typedef Basis<BasisFunctionType> Basis;

    const int trialElementCount = trialElementIndices.size();
    // const int componentCount = m_integral->resultDimension();

    result.resize(trialElementCount);

    // Select integrators
    typedef std::pair<const Integrator*, const Basis*> QuadVariant;
    std::vector<QuadVariant> quadVariants(trialElementCount);
    for (int i = 0; i < trialElementCount; ++i) {
        const int activeTrialElementIndex = trialElementIndices[i];
        const Integrator* integrator =
                &selectIntegrator(pointIndex, activeTrialElementIndex,
                                  nominalDistance);
        quadVariants[i] = QuadVariant(integrator,
                                      (*m_trialBases)[activeTrialElementIndex]);
    }

    // Integration will proceed in batches of test elements having the same
    // "quadrature variant", i.e. integrator and basis

    // Find all the unique quadrature variants present
    typedef std::set<QuadVariant> QuadVariantSet;
    // Set of unique quadrature variants
    QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

    std::vector<int> activeTrialElementIndices;
    activeTrialElementIndices.reserve(trialElementCount);
    std::vector<arma::Mat<ResultType>*> activeLocalResults;
    activeLocalResults.reserve(trialElementCount);

    // Now loop over unique quadrature variants
    for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
         it != uniqueQuadVariants.end(); ++it) {
        const QuadVariant activeQuadVariant = *it;
        const Integrator& activeIntegrator = *activeQuadVariant.first;
        const Basis& activeTrialBasis = *activeQuadVariant.second;

        // Find all the test elements for which quadrature should proceed
        // according to the current quadrature variant
        activeTrialElementIndices.clear();
        activeLocalResults.clear();
        for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
            if (quadVariants[trialIndex] == activeQuadVariant) {
                activeTrialElementIndices.push_back(trialElementIndices[trialIndex]);
                activeLocalResults.push_back(&result[trialIndex]);
            }

        // Integrate!
        activeIntegrator.integrate(pointIndex, componentIndex,
                                   activeTrialElementIndices,
                                   activeTrialBasis,
                                   activeLocalResults);
    }
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
evaluateLocalContributions(
        const std::vector<int>& pointIndices,
        const std::vector<int>& trialElementIndices,
        Fiber::_2dArray<arma::Mat<ResultType> >& result,
        CoordinateType nominalDistance)
{
    typedef Fiber::Basis<BasisFunctionType> Basis;

    const int pointCount = pointIndices.size();
    const int trialElementCount = trialElementIndices.size();
    result.set_size(pointCount, trialElementCount);

    // Find cached matrices; select integrators to calculate non-cached ones
    typedef std::pair<const Integrator*, const Basis*> QuadVariant;
    Fiber::_2dArray<QuadVariant> quadVariants(pointCount, trialElementCount);

    for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
        for (int pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
            const int activePointIndex = pointIndices[pointIndex];
            const int activeTrialElementIndex = trialElementIndices[trialIndex];
            const Integrator* integrator =
                    &selectIntegrator(activePointIndex,
                                      activeTrialElementIndex, nominalDistance);
            quadVariants(pointIndex, trialIndex) = QuadVariant(
                        integrator, (*m_trialBases)[activeTrialElementIndex]);
        }

    // Integration will proceed in batches of element pairs having the same
    // "quadrature variant", i.e. integrator, test basis and trial basis

    // Find all the unique quadrature variants present
    typedef std::set<QuadVariant> QuadVariantSet;
    // Set of unique quadrature variants
    QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

    typedef typename KernelTrialIntegrator<
            BasisFunctionType, KernelType, ResultType>::PointElementIndexPair
            PointElementIndexPair;
    std::vector<PointElementIndexPair> activePointElementPairs;
    std::vector<arma::Mat<ResultType>*> activeLocalResults;
    activePointElementPairs.reserve(pointCount * trialElementCount);
    activeLocalResults.reserve(pointCount * trialElementCount);

    // Now loop over unique quadrature variants
    for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
         it != uniqueQuadVariants.end(); ++it) {
        const QuadVariant activeQuadVariant = *it;
        const Integrator& activeIntegrator = *it->first;
        const Basis& activeTrialBasis = *it->second;

        // Find all the element pairs for which quadrature should proceed
        // according to the current quadrature variant
        activePointElementPairs.clear();
        activeLocalResults.clear();
        for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
            for (int pointIndex = 0; pointIndex < pointCount; ++pointIndex)
            if (quadVariants(pointIndex, trialIndex) == activeQuadVariant) {
                    activePointElementPairs.push_back(
                                PointElementIndexPair(
                                    pointIndices[pointIndex],
                                    trialElementIndices[trialIndex]));
                    activeLocalResults.push_back(&result(pointIndex, trialIndex));
                }

        // Integrate!
        activeIntegrator.integrate(activePointElementPairs,
                                   activeTrialBasis, activeLocalResults);
    }
}


template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
int
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::resultDimension() const
{
    return m_integral->resultDimension();
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
precalculateElementSizesAndCenters()
{
    precalculateElementSizesAndCentersForSingleGrid(
                *m_rawGeometry,
                m_elementSizesSquared, m_elementCenters,
                m_averageElementSize);
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
void
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
precalculateElementSizesAndCentersForSingleGrid(
        const RawGridGeometry<CoordinateType>& rawGeometry,
        std::vector<CoordinateType>& elementSizesSquared,
        arma::Mat<CoordinateType>& elementCenters,
        CoordinateType& averageElementSize) const
{
    const size_t elementCount = rawGeometry.elementCount();
    const int worldDim = rawGeometry.worldDimension();

    averageElementSize = 0.; // We will store here temporarily
                             // the sum of element sizes
    elementSizesSquared.resize(elementCount);
    for (int e = 0; e < elementCount; ++e) {
        elementSizesSquared[e] = elementSizeSquared(e, rawGeometry);
        averageElementSize += sqrt(elementSizesSquared[e]);
    }
    averageElementSize /= elementCount;

    elementCenters.set_size(worldDim, elementCount);
    for (int e = 0; e < elementCount; ++e)
        elementCenters.col(e) = elementCenter(e, rawGeometry);
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
const KernelTrialIntegrator<BasisFunctionType, KernelType, ResultType>&
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
selectIntegrator(int pointIndex, int trialElementIndex,
                 CoordinateType nominalDistance)
{
    SingleQuadratureDescriptor desc;
    desc.vertexCount = m_rawGeometry->elementCornerCount(trialElementIndex);
    desc.order = order(pointIndex, trialElementIndex, nominalDistance);
    return getIntegrator(desc);
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
int
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::order(
        int pointIndex, int trialElementIndex,
        CoordinateType nominalDistance) const
{
    // Order required for exact quadrature on affine elements with a constant kernel
    int trialBasisOrder = (*m_trialBases)[trialElementIndex]->order();
    int defaultQuadOrder = 2 * trialBasisOrder;

    CoordinateType normalisedDistance;
    if (nominalDistance < 0.) {
        CoordinateType elementSizeSquared =
                m_elementSizesSquared[trialElementIndex];
        CoordinateType distanceSquared =
                pointElementDistanceSquared(pointIndex, trialElementIndex);
        CoordinateType normalisedDistanceSquared =
                distanceSquared / elementSizeSquared;
        normalisedDistance = sqrt(normalisedDistanceSquared);
    } else
        normalisedDistance = nominalDistance / m_averageElementSize;

    const QuadratureOptions& options =
            m_accuracyOptions.singleRegular(); // (normalisedDistance); // TODO
    return options.quadratureOrder(defaultQuadOrder);
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
inline
typename DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
BasisFunctionType, KernelType, ResultType, GeometryFactory>::CoordinateType
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
BasisFunctionType, KernelType, ResultType, GeometryFactory>::elementSizeSquared(
        int elementIndex, const RawGridGeometry<CoordinateType>& rawGeometry) const
{
    // This implementation could be optimised
    CoordinateType maxEdgeLengthSquared = 0.;
    const arma::Mat<int>& cornerIndices = rawGeometry.elementCornerIndices();
    const arma::Mat<CoordinateType>& vertices = rawGeometry.vertices();
    arma::Col<CoordinateType> edge;
    if (cornerIndices(cornerIndices.n_rows - 1, elementIndex) == -1) {
        // Triangular element
        const int cornerCount = 3;
        for (int i = 0; i < cornerCount; ++i) {
            edge = vertices.col(cornerIndices((i + 1) % cornerCount, elementIndex)) -
                    vertices.col(cornerIndices(i, elementIndex));
            CoordinateType edgeLengthSquared = arma::dot(edge, edge);
            maxEdgeLengthSquared = std::max(maxEdgeLengthSquared, edgeLengthSquared);
        }
    } else {
        // Quadrilateral element. We assume it is convex.
        edge = vertices.col(cornerIndices(2, elementIndex)) -
                vertices.col(cornerIndices(0, elementIndex));
        maxEdgeLengthSquared = arma::dot(edge, edge);
        edge = vertices.col(cornerIndices(3, elementIndex)) -
                vertices.col(cornerIndices(1, elementIndex));
        CoordinateType edgeLengthSquared = arma::dot(edge, edge);
        maxEdgeLengthSquared = std::max(maxEdgeLengthSquared, edgeLengthSquared);
    }
    return maxEdgeLengthSquared;
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
inline
arma::Col<typename DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
BasisFunctionType, KernelType, ResultType, GeometryFactory>::CoordinateType>
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
BasisFunctionType, KernelType, ResultType, GeometryFactory>::elementCenter(
        int elementIndex, const RawGridGeometry<CoordinateType>& rawGeometry) const
{
    const arma::Mat<int>& cornerIndices = rawGeometry.elementCornerIndices();
    const arma::Mat<CoordinateType>& vertices = rawGeometry.vertices();
    const int maxCornerCount = cornerIndices.n_rows;
    // each element has at least one corner
    arma::Col<CoordinateType> center(vertices.col(cornerIndices(0, elementIndex)));
    int i = 1;
    for (; i < maxCornerCount; ++i) {
        int cornerIndex = cornerIndices(i, elementIndex);
        if (cornerIndex == -1)
            break;
        center += vertices.col(cornerIndex);
    }
    // now i contains the number of corners of the specified element
    center /= i;
    return center;
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
inline
typename DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
BasisFunctionType, KernelType, ResultType, GeometryFactory>::CoordinateType
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
BasisFunctionType, KernelType, ResultType, GeometryFactory>::pointElementDistanceSquared(
        int pointIndex, int trialElementIndex) const
{
    arma::Col<CoordinateType> diff =
            m_points.col(pointIndex) -
            m_elementCenters.col(trialElementIndex);
    return arma::dot(diff, diff);
}

template <typename BasisFunctionType, typename KernelType,
          typename ResultType, typename GeometryFactory>
const KernelTrialIntegrator<BasisFunctionType, KernelType, ResultType>&
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<BasisFunctionType,
KernelType, ResultType, GeometryFactory>::
getIntegrator(const SingleQuadratureDescriptor& desc)
{
    typename IntegratorMap::const_iterator it = m_kernelTrialIntegrators.find(desc);
    // Note: as far as I understand TBB's docs, .end() keeps pointing to the
    // same element even if another thread inserts a new element into the map
    if (it != m_kernelTrialIntegrators.end()) {
        //std::cout << "getIntegrator(: " << desc << "): integrator found" << std::endl;
        return *it->second;
    }
    //std::cout << "getIntegrator(: " << desc << "): integrator not found" << std::endl;

    // Integrator doesn't exist yet and must be created.
    Integrator* integrator = 0;
    // Create a quadrature rule
    arma::Mat<CoordinateType> trialPoints;
    std::vector<CoordinateType> trialWeights;

    fillSingleQuadraturePointsAndWeights(desc.vertexCount,
                                         desc.order,
                                         trialPoints, trialWeights);
    typedef NumericalKernelTrialIntegrator<BasisFunctionType,
            KernelType, ResultType, GeometryFactory> ConcreteIntegrator;
    integrator = new ConcreteIntegrator(
                trialPoints, trialWeights,
                m_points,
                *m_geometryFactory, *m_rawGeometry,
                *m_kernels, *m_trialTransformations, *m_integral);

    // Attempt to insert the newly created integrator into the map
    std::pair<typename IntegratorMap::iterator, bool> result =
            m_kernelTrialIntegrators.insert(std::make_pair(desc, integrator));
    if (result.second)
        // Insertion succeeded. The newly created integrator will be deleted in
        // our own destructor
        ;
    else
        // Insertion failed -- another thread was faster. Delete the newly
        // created integrator.
        delete integrator;

    //    if (result.second)
    //        std::cout << "getIntegrator(: " << desc << "): insertion succeeded" << std::endl;
    //    else
    //        std::cout << "getIntegrator(: " << desc << "): insertion failed" << std::endl;

    // Return pointer to the integrator that ended up in the map.
    return *result.first->second;
}

} // namespace Fiber
