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

#ifndef fiber_default_local_assembler_for_integral_operators_on_surfaces_hpp
#define fiber_default_local_assembler_for_integral_operators_on_surfaces_hpp

#include "../common/common.hpp"

#include "local_assembler_for_integral_operators.hpp"

#include "_2d_array.hpp"
#include "accuracy_options.hpp"
#include "default_local_assembler_for_operators_on_surfaces_utilities.hpp"
#include "element_pair_topology.hpp"
#include "numerical_quadrature.hpp"
#include "parallelization_options.hpp"
#include "shared_ptr.hpp"
#include "test_kernel_trial_integrator.hpp"
#include "verbosity_level.hpp"

#include <boost/static_assert.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/mutex.h>
#include <cstring>
#include <climits>
#include <set>
#include <utility>
#include <vector>

namespace Fiber {

/** \cond FORWARD_DECL */
class OpenClHandler;
template <typename CoordinateType> class CollectionOfShapesetTransformations;
template <typename ValueType> class CollectionOfKernels;
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class TestKernelTrialIntegral;

template <typename CoordinateType> class RawGridGeometry;

template <typename CoordinateType>
class QuadratureDescriptorSelectorForIntegralOperators;
template <typename CoordinateType> class DoubleQuadratureRuleFamily;
/** \endcond */

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
class DefaultLocalAssemblerForIntegralOperatorsOnSurfaces
    : public LocalAssemblerForIntegralOperators<ResultType> {
public:
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  DefaultLocalAssemblerForIntegralOperatorsOnSurfaces(
      const shared_ptr<const GeometryFactory> &testGeometryFactory,
      const shared_ptr<const GeometryFactory> &trialGeometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &testRawGeometry,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &trialRawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const CollectionOfKernels<KernelType>> &kernel,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const TestKernelTrialIntegral<
          BasisFunctionType, KernelType, ResultType>> &integral,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions,
      VerbosityLevel::Level verbosityLevel, bool cacheSingularIntegrals,
      const shared_ptr<const QuadratureDescriptorSelectorForIntegralOperators<
          CoordinateType>> &quadDescSelector,
      const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType>> &
          quadRuleFamily);
  virtual ~DefaultLocalAssemblerForIntegralOperatorsOnSurfaces();

public:
  virtual void
  evaluateLocalWeakForms(CallVariant callVariant,
                         const std::vector<int> &elementIndicesA,
                         int elementIndexB, LocalDofIndex localDofIndexB,
                         std::vector<arma::Mat<ResultType>> &result,
                         CoordinateType nominalDistance = -1.);

  virtual void
  evaluateLocalWeakForms(const std::vector<int> &testElementIndices,
                         const std::vector<int> &trialElementIndices,
                         Fiber::_2dArray<arma::Mat<ResultType>> &result,
                         CoordinateType nominalDistance = -1.);

  virtual CoordinateType estimateRelativeScale(CoordinateType minDist) const;

private:
  /** \cond PRIVATE */
  typedef TestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType>
  Integrator;
  typedef typename Integrator::ElementIndexPair ElementIndexPair;
  typedef DefaultLocalAssemblerForOperatorsOnSurfacesUtilities<
      BasisFunctionType> Utilities;

  /** \brief Alternative comparison functor for pairs.
   *
   *  This functor can be used to sort pairs according to the second member
   *  and then, in case of equality, according to the first member. */
  template <typename T1, typename T2> struct alternative_less {
    bool operator()(const std::pair<T1, T2> &a,
                    const std::pair<T1, T2> &b) const {
      return a.second < b.second ||
             (!(b.second < a.second) && a.first < b.first);
    }
  };

  /** \brief Comparison functor for element-index pairs.
   *
   *  This functor sorts element index pairs first after the trial element
   *  index (second member) and then, in case of equality, after the test
   *  element index (first member) */
  typedef alternative_less<typename ElementIndexPair::first_type,
                           typename ElementIndexPair::second_type>
  ElementIndexPairCompare;

  /** \brief Set of element index pairs.
   *
   *  The alternative sorting (first after the trial element index) is used
   *  because profiling has shown that evaluateLocalWeakForms is called more
   *  often in the TEST_TRIAL mode (with a single trial element index) than
   *  in the TRIAL_TEST mode. Therefore the singular integral cache is
   *  indexed with trial element index, and this sorting mode makes it easier
   *  to construct such cache. */
  typedef std::set<ElementIndexPair, ElementIndexPairCompare>
  ElementIndexPairSet;

  bool testAndTrialGridsAreIdentical() const;

  void cacheSingularLocalWeakForms();
  void findPairsOfAdjacentElements(ElementIndexPairSet &pairs) const;
  void cacheLocalWeakForms(const ElementIndexPairSet &elementIndexPairs);

  const Integrator &selectIntegrator(int testElementIndex,
                                     int trialElementIndex,
                                     CoordinateType nominalDistance = -1.);

  const Integrator &getIntegrator(const DoubleQuadratureDescriptor &index);

private:
  shared_ptr<const GeometryFactory> m_testGeometryFactory;
  shared_ptr<const GeometryFactory> m_trialGeometryFactory;
  shared_ptr<const RawGridGeometry<CoordinateType>> m_testRawGeometry;
  shared_ptr<const RawGridGeometry<CoordinateType>> m_trialRawGeometry;
  shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
  m_testShapesets;
  shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
  m_trialShapesets;
  shared_ptr<const CollectionOfShapesetTransformations<CoordinateType>>
  m_testTransformations;
  shared_ptr<const CollectionOfKernels<KernelType>> m_kernels;
  shared_ptr<const CollectionOfShapesetTransformations<CoordinateType>>
  m_trialTransformations;
  shared_ptr<const TestKernelTrialIntegral<BasisFunctionType, KernelType,
                                           ResultType>> m_integral;
  shared_ptr<const OpenClHandler> m_openClHandler;
  ParallelizationOptions m_parallelizationOptions;
  VerbosityLevel::Level m_verbosityLevel;
  shared_ptr<const QuadratureDescriptorSelectorForIntegralOperators<
      CoordinateType>> m_quadDescSelector;
  shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType>> m_quadRuleFamily;

  typedef tbb::concurrent_unordered_map<DoubleQuadratureDescriptor,
                                        Integrator *> IntegratorMap;
  IntegratorMap m_testKernelTrialIntegrators;
  mutable tbb::mutex m_integratorCreationMutex;

  enum {
    INVALID_INDEX = INT_MAX
  };
  typedef _2dArray<std::pair<int, arma::Mat<ResultType>>> Cache;
  /** \brief Singular integral cache.
   *
   *  This cache stores the preevaluated local weak forms expressed by
   *  singular integrals. A particular item it stored in r'th row and c'th
   *  column stores, in its second member, the local weak form calculated for
   *  the test element with index it.first and the trial element with index
   *  c. In each column, the items are sorted after increasing test element
   *  index. At the end of each column there can be unused items with test
   *  element index set to INVALID_INDEX (= INT_MAX, so that the sorting is
   *  preserved). */
  Cache m_cache;
  /** \endcond */
};

} // namespace Fiber

#include "default_local_assembler_for_integral_operators_on_surfaces_imp.hpp"

#endif
