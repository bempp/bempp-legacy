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

#include "cuda_local_assembler_for_integral_operators_on_surfaces.hpp"

namespace Fiber {

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
CudaLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType, GeometryFactory>::
    CudaLocalAssemblerForIntegralOperatorsOnSurfaces() { }

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
CudaLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType,
    GeometryFactory>::~CudaLocalAssemblerForIntegralOperatorsOnSurfaces() { }

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
void CudaLocalAssemblerForIntegralOperatorsOnSurfaces<
    BasisFunctionType, KernelType, ResultType, GeometryFactory>::
    evaluateLocalWeakForms(CallVariant callVariant,
                           const std::vector<int> &elementIndicesA,
                           int elementIndexB, LocalDofIndex localDofIndexB,
                           std::vector<Matrix<ResultType>> &result,
                           CoordinateType nominalDistance) {

//  typedef Shapeset<BasisFunctionType> Shapeset;
//
//  const int elementACount = elementIndicesA.size();
//  result.resize(elementACount);
//
//  // Get shapesets
//  const std::vector<const Shapeset *> &m_basesA =
//      callVariant == TEST_TRIAL ? *m_testShapesets : *m_trialShapesets;
//  const std::vector<const Shapeset *> &m_basesB =
//      callVariant == TEST_TRIAL ? *m_trialShapesets : *m_testShapesets;
//  std::vector<const Shapeset *> basesA(elementACount);
//  for (int i = 0; i < elementACount; ++i)
//    basesA[i] = m_basesA[elementIndicesA[i]];
//  const Shapeset &basisB = *m_basesB[elementIndexB];
//
//  // Find cached matrices; select integrators to calculate non-cached ones
//  typedef std::pair<const Integrator *, const Shapeset *> QuadVariant;
//  const QuadVariant CACHED(0, 0);
//  std::vector<QuadVariant> quadVariants(elementACount);
//  for (int i = 0; i < elementACount; ++i) {
//    // Try to find matrix in cache
//    const Matrix<ResultType> *cachedLocalWeakForm = 0;
//    if (callVariant == TEST_TRIAL) {
//      const int testElementIndex = elementIndicesA[i];
//      const int trialElementIndex = elementIndexB;
//      for (size_t n = 0; n < m_cache.extent(0); ++n)
//        if (m_cache(n, trialElementIndex).first == testElementIndex) {
//          cachedLocalWeakForm = &m_cache(n, trialElementIndex).second;
//          break;
//        }
//    } else {
//      const int testElementIndex = elementIndexB;
//      const int trialElementIndex = elementIndicesA[i];
//      for (size_t n = 0; n < m_cache.extent(0); ++n)
//        if (m_cache(n, trialElementIndex).first == testElementIndex) {
//          cachedLocalWeakForm = &m_cache(n, trialElementIndex).second;
//          break;
//        }
//    }
//
//    if (cachedLocalWeakForm) { // Matrix found in cache
//      quadVariants[i] = CACHED;
//      if (localDofIndexB == ALL_DOFS)
//        result[i] = *cachedLocalWeakForm;
//      else {
//        if (callVariant == TEST_TRIAL)
//          result[i] = cachedLocalWeakForm->col(localDofIndexB);
//        else
//          result[i] = cachedLocalWeakForm->row(localDofIndexB);
//      }
//    } else {
//      const Integrator *integrator =
//          callVariant == TEST_TRIAL
//              ? &selectIntegrator(elementIndicesA[i], elementIndexB,
//                                  nominalDistance)
//              : &selectIntegrator(elementIndexB, elementIndicesA[i],
//                                  nominalDistance);
//      quadVariants[i] = QuadVariant(integrator, basesA[i]);
//    }
//  }
//
//  // Integration will proceed in batches of test elements having the same
//  // "quadrature variant", i.e. integrator and shapeset
//
//  // Find all the unique quadrature variants present
//  typedef std::set<QuadVariant> QuadVariantSet;
//  // Set of unique quadrature variants
//  QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());
//
//  std::vector<int> activeElementIndicesA;
//  activeElementIndicesA.reserve(elementACount);
//  std::vector<Matrix<ResultType> *> activeLocalResults;
//  activeLocalResults.reserve(elementACount);
//
//  // Now loop over unique quadrature variants
//  for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
//       it != uniqueQuadVariants.end(); ++it) {
//    const QuadVariant activeQuadVariant = *it;
//    if (activeQuadVariant == CACHED)
//      continue;
//    const Integrator &activeIntegrator = *it->first;
//    const Shapeset &activeBasisA = *it->second;
//
//    // Find all the test elements for which quadrature should proceed
//    // according to the current quadrature variant
//    activeElementIndicesA.clear();
//    activeLocalResults.clear();
//    for (int indexA = 0; indexA < elementACount; ++indexA)
//      if (quadVariants[indexA] == activeQuadVariant) {
//        activeElementIndicesA.push_back(elementIndicesA[indexA]);
//        activeLocalResults.push_back(&result[indexA]);
//      }
//
//    // Integrate!
//    activeIntegrator.integrate(callVariant, activeElementIndicesA,
//                               elementIndexB, activeBasisA, basisB,
//                               localDofIndexB, activeLocalResults);
//  }
}

} // namespace Fiber
