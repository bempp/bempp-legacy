// Copyright (C) 2011-2014 by the BEM++ Authors
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

#ifndef bempp_weak_form_hmat_assembly_helper_hpp
#define bempp_weak_form_hmat_assembly_helper_hpp

#include "../common/common.hpp"

#include "../common/armadillo_fwd.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/types.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../hmat/common.hpp"

#include <tbb/atomic.h>
#include <tbb/concurrent_unordered_map.h>
#include <vector>

namespace hmat {
/** \cond FORWARD_DECL */
template <int N = 2> class BlockClusterTree;
template <int N = 2> class BlockClusterTreeNode;

/** \endcond */
}

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForIntegralOperators;
/** \endcond */

} // namespace Fiber

namespace Bempp {

/** \cond FORWARD_DECL */
class AssemblyOptions;
template <typename ResultType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType> class LocalDofListsCache;
template <typename BasisFunctionType> class Space;
/** \endcond */

/** \ingroup weak_form_assembly_internal
 *  \brief Class whose methods are called by Ahmed during assembly in the ACA
 * mode.
 */
template <typename BasisFunctionType, typename ResultType>
class WeakFormHMatAssemblyHelper {
public:
  typedef DiscreteBoundaryOperator<ResultType> DiscreteLinOp;
  typedef Fiber::LocalAssemblerForIntegralOperators<ResultType> LocalAssembler;
  typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;
  typedef typename Fiber::ScalarTraits<ResultType>::RealType MagnitudeType;

  WeakFormHMatAssemblyHelper(
      const Space<BasisFunctionType> &testSpace,
      const Space<BasisFunctionType> &trialSpace,
      const shared_ptr<hmat::BlockClusterTree<2>> blockClusterTree,
      const std::vector<LocalAssembler *> &assemblers,
      const std::vector<const DiscreteLinOp *> &sparseTermsToAdd,
      const std::vector<ResultType> &denseTermsMultipliers,
      const std::vector<ResultType> &sparseTermsMultipliers,
      const AssemblyOptions &options);

  /** \brief Evaluate entries of a general block. */

  //  void
  //  computeMatrixBlock(const hmat::IndexRangeType &testIndexRange,
  //                     const hmat::IndexRangeType &trialIndexRange,
  //                     arma::Mat<ResultType> &data,
  //                     const hmat::BlockClusterTreeNode &blockClusterTreeNode,
  //                     bool countAccessedEntries = true) const;
  //
  //  /** \brief Return the number of entries in the matrix that have been
  //   *  accessed so far. */
  //  size_t accessedEntryCount() const;
  //
  //  /** \brief Reset the number of entries in the matrix that have been
  //   *  accessed so far. */
  //  void resetAccessedEntryCount();

private:
  //  MagnitudeType estimateMinimumDistance(
  //      const hmat::BlockClusterTreeNode &blockClusterTreeNode) const;

private:
  /** \cond PRIVATE */
  const Space<BasisFunctionType> &m_testSpace;
  const Space<BasisFunctionType> &m_trialSpace;
  const shared_ptr<const hmat::BlockClusterTree<>> m_blockClusterTree;
  const std::vector<LocalAssembler *> &m_assemblers;
  const std::vector<const DiscreteLinOp *> &m_sparseTermsToAdd;
  const std::vector<ResultType> &m_denseTermsMultipliers;
  const std::vector<ResultType> &m_sparseTermsMultipliers;
  const AssemblyOptions &m_options;
  bool m_indexWithGlobalDofs;
  bool m_uniformQuadratureOrder;

  shared_ptr<LocalDofListsCache<BasisFunctionType>> m_testDofListsCache,
      m_trialDofListsCache;

  mutable tbb::atomic<size_t> m_accessedEntryCount;

  typedef tbb::concurrent_unordered_map<
      shared_ptr<const hmat::BlockClusterTreeNode<>>, CoordinateType,
      std::hash<shared_ptr<const hmat::BlockClusterTreeNode<>>>> DistanceMap;
  mutable DistanceMap m_distancesCache;

  /** \endcond */
};

} // namespace Bempp

#endif
