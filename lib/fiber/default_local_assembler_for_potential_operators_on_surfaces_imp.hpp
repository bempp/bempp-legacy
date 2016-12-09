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

#include "kernel_trial_integral.hpp"
#include "numerical_kernel_trial_integrator.hpp"
#include "quadrature_descriptor_selector_for_potential_operators.hpp"
#include "raw_grid_geometry.hpp"
#include "serial_blas_region.hpp"
#include "shapeset.hpp"
#include "single_quadrature_rule_family.hpp"

#include <cassert>
#include <tbb/parallel_for.h>

#include "../common/auto_timer.hpp"

namespace Fiber {

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType, GeometryFactory>::
    DefaultLocalAssemblerForPotentialOperatorsOnSurfaces(
        const Matrix<CoordinateType> &points,
        const shared_ptr<const GeometryFactory> &geometryFactory,
        const shared_ptr<const RawGridGeometry<CoordinateType>> &rawGeometry,
        const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
            &trialShapesets,
        const shared_ptr<const CollectionOfKernels<KernelType>> &kernels,
        const shared_ptr<const CollectionOfShapesetTransformations<
            CoordinateType>> &trialTransformations,
        const shared_ptr<const KernelTrialIntegral<
            BasisFunctionType, KernelType, ResultType>> &integral,
        const ParallelizationOptions &parallelizationOptions,
        VerbosityLevel::Level verbosityLevel,
        const shared_ptr<
            const QuadratureDescriptorSelectorForPotentialOperators<
                BasisFunctionType>> &quadDescSelector,
        const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>>
            &quadRuleFamily)
    : m_points(points), m_geometryFactory(geometryFactory),
      m_rawGeometry(rawGeometry), m_trialShapesets(trialShapesets),
      m_kernels(kernels), m_trialTransformations(trialTransformations),
      m_integral(integral), m_parallelizationOptions(parallelizationOptions),
      m_verbosityLevel(verbosityLevel), m_quadDescSelector(quadDescSelector),
      m_quadRuleFamily(quadRuleFamily) {
  Utilities::checkConsistencyOfGeometryAndShapesets(*rawGeometry,
                                                    *trialShapesets);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::~DefaultLocalAssemblerForPotentialOperatorsOnSurfaces() {
  // Note: obviously the destructor is assumed to be called only after
  // all threads have ceased using the assembler!

  for (typename IntegratorMap::const_iterator it =
           m_kernelTrialIntegrators.begin();
       it != m_kernelTrialIntegrators.end(); ++it)
    delete it->second;
  m_kernelTrialIntegrators.clear();
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType, GeometryFactory>::
    evaluateLocalContributions(const std::vector<int> &pointIndices,
                               int trialElementIndex,
                               LocalDofIndex localTrialDofIndex,
                               std::vector<Matrix<ResultType>> &result,
                               CoordinateType nominalDistance) {
  typedef Shapeset<BasisFunctionType> Shapeset;

  const int pointCount = pointIndices.size();
  // const int componentCount = m_integral->resultDimension();

  // Get shapeset
  const Shapeset &trialShapeset = *((*m_trialShapesets)[trialElementIndex]);
  assert(
      localTrialDofIndex == ALL_DOFS ||
      (localTrialDofIndex >= 0 && localTrialDofIndex < trialShapeset.size()));
  result.resize(pointCount);

  // Select integrators
  typedef const Integrator *QuadVariant;
  std::vector<QuadVariant> quadVariants(pointCount);
  for (int i = 0; i < pointCount; ++i) {
    const Integrator *integrator =
        &selectIntegrator(pointIndices[i], trialElementIndex, nominalDistance);
    quadVariants[i] = QuadVariant(integrator);
  }

  // Integration will proceed in batches of test elements having the same
  // "quadrature variant", i.e. integrator and shapeset

  // Find all the unique quadrature variants present
  typedef std::set<QuadVariant> QuadVariantSet;
  // Set of unique quadrature variants
  QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

  std::vector<int> activePointIndices;
  activePointIndices.reserve(pointCount);
  std::vector<Matrix<ResultType> *> activeLocalResults;
  activeLocalResults.reserve(pointCount);

  // Now loop over unique quadrature variants
  for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
       it != uniqueQuadVariants.end(); ++it) {
    const QuadVariant activeQuadVariant = *it;
    const Integrator &activeIntegrator = *activeQuadVariant;

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
                               trialShapeset, localTrialDofIndex,
                               activeLocalResults);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType, GeometryFactory>::
    evaluateLocalContributions(int pointIndex, int componentIndex,
                               const std::vector<int> &trialElementIndices,
                               std::vector<Matrix<ResultType>> &result,
                               CoordinateType nominalDistance) {
  typedef Shapeset<BasisFunctionType> Shapeset;

  const int trialElementCount = trialElementIndices.size();
  // const int componentCount = m_integral->resultDimension();

  result.resize(trialElementCount);

  // Select integrators
  typedef std::pair<const Integrator *, const Shapeset *> QuadVariant;
  std::vector<QuadVariant> quadVariants(trialElementCount);
  for (int i = 0; i < trialElementCount; ++i) {
    const int activeTrialElementIndex = trialElementIndices[i];
    const Integrator *integrator =
        &selectIntegrator(pointIndex, activeTrialElementIndex, nominalDistance);
    quadVariants[i] =
        QuadVariant(integrator, (*m_trialShapesets)[activeTrialElementIndex]);
  }

  // Integration will proceed in batches of test elements having the same
  // "quadrature variant", i.e. integrator and shapeset

  // Find all the unique quadrature variants present
  typedef std::set<QuadVariant> QuadVariantSet;
  // Set of unique quadrature variants
  QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

  std::vector<int> activeTrialElementIndices;
  activeTrialElementIndices.reserve(trialElementCount);
  std::vector<Matrix<ResultType> *> activeLocalResults;
  activeLocalResults.reserve(trialElementCount);

  // Now loop over unique quadrature variants
  for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
       it != uniqueQuadVariants.end(); ++it) {
    const QuadVariant activeQuadVariant = *it;
    const Integrator &activeIntegrator = *activeQuadVariant.first;
    const Shapeset &activeTrialShapeset = *activeQuadVariant.second;

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
                               activeTrialElementIndices, activeTrialShapeset,
                               activeLocalResults);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType, GeometryFactory>::
    evaluateLocalContributions(const std::vector<int> &pointIndices,
                               const std::vector<int> &trialElementIndices,
                               Fiber::_2dArray<Matrix<ResultType>> &result,
                               CoordinateType nominalDistance) {
  typedef Fiber::Shapeset<BasisFunctionType> Shapeset;

  const int pointCount = pointIndices.size();
  const int trialElementCount = trialElementIndices.size();
  result.set_size(pointCount, trialElementCount);

  // Find cached matrices; select integrators to calculate non-cached ones
  typedef std::pair<const Integrator *, const Shapeset *> QuadVariant;
  Fiber::_2dArray<QuadVariant> quadVariants(pointCount, trialElementCount);

  for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
    for (int pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
      const int activePointIndex = pointIndices[pointIndex];
      const int activeTrialElementIndex = trialElementIndices[trialIndex];
      const Integrator *integrator = &selectIntegrator(
          activePointIndex, activeTrialElementIndex, nominalDistance);
      quadVariants(pointIndex, trialIndex) =
          QuadVariant(integrator, (*m_trialShapesets)[activeTrialElementIndex]);
    }

  // Integration will proceed in batches of element pairs having the same
  // "quadrature variant", i.e. integrator, test shapeset and trial shapeset

  // Find all the unique quadrature variants present
  typedef std::set<QuadVariant> QuadVariantSet;
  // Set of unique quadrature variants
  QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

  typedef typename KernelTrialIntegrator<BasisFunctionType, KernelType,
                                         ResultType>::PointElementIndexPair
      PointElementIndexPair;
  std::vector<PointElementIndexPair> activePointElementPairs;
  std::vector<Matrix<ResultType> *> activeLocalResults;
  activePointElementPairs.reserve(pointCount * trialElementCount);
  activeLocalResults.reserve(pointCount * trialElementCount);

  // Now loop over unique quadrature variants
  for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
       it != uniqueQuadVariants.end(); ++it) {
    const QuadVariant activeQuadVariant = *it;
    const Integrator &activeIntegrator = *it->first;
    const Shapeset &activeTrialShapeset = *it->second;

    // Find all the element pairs for which quadrature should proceed
    // according to the current quadrature variant
    activePointElementPairs.clear();
    activeLocalResults.clear();
    for (int trialIndex = 0; trialIndex < trialElementCount; ++trialIndex)
      for (int pointIndex = 0; pointIndex < pointCount; ++pointIndex)
        if (quadVariants(pointIndex, trialIndex) == activeQuadVariant) {
          activePointElementPairs.push_back(PointElementIndexPair(
              pointIndices[pointIndex], trialElementIndices[trialIndex]));
          activeLocalResults.push_back(&result(pointIndex, trialIndex));
        }

    // Integrate!
    activeIntegrator.integrate(activePointElementPairs, activeTrialShapeset,
                               activeLocalResults);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
int DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::resultDimension() const {
  return m_integral->resultDimension();
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
typename DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType, GeometryFactory>::CoordinateType
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::estimateRelativeScale(CoordinateType minDist) const {
  return m_kernels->estimateRelativeScale(minDist);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
const KernelTrialIntegrator<BasisFunctionType, KernelType, ResultType> &
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::selectIntegrator(int pointIndex, int trialElementIndex,
                                       CoordinateType nominalDistance) {
  const Eigen::Map<Vector<CoordinateType>> pointCoords(
      m_points.col(pointIndex).data(), m_points.rows());
  // m_points.unsafe_col(pointIndex);
  SingleQuadratureDescriptor desc = m_quadDescSelector->quadratureDescriptor(
      pointCoords, trialElementIndex, nominalDistance);
  return getIntegrator(desc);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
const KernelTrialIntegrator<BasisFunctionType, KernelType, ResultType> &
DefaultLocalAssemblerForPotentialOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::getIntegrator(const SingleQuadratureDescriptor &desc) {
  typename IntegratorMap::const_iterator it =
      m_kernelTrialIntegrators.find(desc);
  // Note: as far as I understand TBB's docs, .end() keeps pointing to the
  // same element even if another thread inserts a new element into the map
  if (it != m_kernelTrialIntegrators.end()) {
    // std::cout << "getIntegrator(: " << desc << "): integrator found" <<
    // std::endl;
    return *it->second;
  }
  // std::cout << "getIntegrator(: " << desc << "): integrator not found" <<
  // std::endl;

  // Integrator doesn't exist yet and must be created.
  Integrator *integrator = 0;
  // Create a quadrature rule
  Matrix<CoordinateType> trialPoints;
  std::vector<CoordinateType> trialWeights;

  m_quadRuleFamily->fillQuadraturePointsAndWeights(desc, trialPoints,
                                                   trialWeights);
  typedef NumericalKernelTrialIntegrator<BasisFunctionType, KernelType,
                                         ResultType,
                                         GeometryFactory> ConcreteIntegrator;
  integrator = new ConcreteIntegrator(
      trialPoints, trialWeights, m_points, *m_geometryFactory, *m_rawGeometry,
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
  //        std::cout << "getIntegrator(: " << desc << "): insertion succeeded"
  // << std::endl;
  //    else
  //        std::cout << "getIntegrator(: " << desc << "): insertion failed" <<
  // std::endl;

  // Return pointer to the integrator that ended up in the map.
  return *result.first->second;
}

} // namespace Fiber
