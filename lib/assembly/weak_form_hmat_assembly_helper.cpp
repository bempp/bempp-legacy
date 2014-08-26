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

#include "weak_form_hmat_assembly_helper.hpp"
#include "assembly_options.hpp"
#include "../hmat/block_cluster_tree.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/types.hpp"
#include "local_dof_lists_cache.hpp"

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
WeakFormHMatAssemblyHelper<BasisFunctionType, ResultType>::
    WeakFormHMatAssemblyHelper(
        const Space<BasisFunctionType> &testSpace,
        const Space<BasisFunctionType> &trialSpace,
        const shared_ptr<hmat::DefaultBlockClusterTreeType> blockClusterTree,
        const std::vector<LocalAssembler *> &assemblers,
        const std::vector<const DiscreteLinOp *> &sparseTermsToAdd,
        const std::vector<ResultType> &denseTermsMultipliers,
        const std::vector<ResultType> &sparseTermsMultipliers,
        const AssemblyOptions &options)
    : m_testSpace(testSpace), m_trialSpace(trialSpace),
      m_blockClusterTree(blockClusterTree), m_assemblers(assemblers),
      m_sparseTermsToAdd(sparseTermsToAdd),
      m_denseTermsMultipliers(denseTermsMultipliers),
      m_sparseTermsMultipliers(sparseTermsMultipliers), m_options(options),
      m_indexWithGlobalDofs(m_options.acaOptions().mode !=
                            AcaOptions::HYBRID_ASSEMBLY),
      m_uniformQuadratureOrder(
          m_options.isQuadratureOrderUniformInEachCluster()),
      m_testDofListsCache(new LocalDofListsCache<BasisFunctionType>(
          m_testSpace,
          blockClusterTree->rowClusterTree()->hMatDofToOriginalDofMap(),
          m_indexWithGlobalDofs)),
      m_trialDofListsCache(new LocalDofListsCache<BasisFunctionType>(
          m_trialSpace,
          blockClusterTree->columnClusterTree()->hMatDofToOriginalDofMap(),
          m_indexWithGlobalDofs)) {

  if (!m_indexWithGlobalDofs && !m_sparseTermsToAdd.empty())
    throw std::invalid_argument(
        "WeakFormAcaAssemblyHelper::WeakFormAcaAssemblyHelper(): "
        "combining sparse and dense terms in hybrid ACA mode "
        "is not supported at present");
  for (size_t i = 0; i < assemblers.size(); ++i)
    if (!assemblers[i])
      throw std::invalid_argument(
          "WeakFormAcaAssemblyHelper::WeakFormAcaAssemblyHelper(): "
          "no elements of the 'assemblers' vector may be null");
  for (size_t i = 0; i < sparseTermsToAdd.size(); ++i)
    if (!sparseTermsToAdd[i])
      throw std::invalid_argument(
          "WeakFormAcaAssemblyHelper::WeakFormAcaAssemblyHelper(): "
          "no elements of the 'sparseTermsToAdd' vector may be null");
  m_accessedEntryCount = 0;
}

template <typename BasisFunctionType, typename ResultType>
typename WeakFormHMatAssemblyHelper<BasisFunctionType,
                                    ResultType>::MagnitudeType
WeakFormHMatAssemblyHelper<BasisFunctionType, ResultType>::
    estimateMinimumDistance(const hmat::DefaultBlockClusterTreeNodeType &
                                blockClusterTreeNode) const {

  return MagnitudeType(
      blockClusterTreeNode.data()
          .rowClusterTreeNode->data()
          .boundingBox.distance(blockClusterTreeNode.data()
                                    .columnClusterTreeNode->data()
                                    .boundingBox));
}

template <typename BasisFunctionType, typename ResultType>
void
WeakFormHMatAssemblyHelper<BasisFunctionType, ResultType>::computeMatrixBlock(
    const hmat::IndexRangeType &testIndexRange,
    const hmat::IndexRangeType &trialIndexRange, arma::Mat<ResultType> &data,
    const hmat::DefaultBlockClusterTreeNodeType &blockClusterTreeNode,
    bool countAccessedEntries) const {

  if (countAccessedEntries)
    m_accessedEntryCount += (testIndexRange[1] - testIndexRange[0]) *
                            (trialIndexRange[1] - trialIndexRange[0]);

	
   const CoordinateType minDist = 
	m_uniformQuadratureOrder ? estimateMinimumDistance(blockClusterTreeNode) : -1;

   shared_ptr<const LocalDofLists<BasisFunctionType>> testDofLists =
       m_testDofListsCache->get(testIndexRange[0],
                               testIndexRange[1] - testIndexRange[0]);
   shared_ptr<const LocalDofLists<BasisFunctionType>> trialDofLists =
       m_trialDofListsCache->get(trialIndexRange[0],
                               trialIndexRange[1] - trialIndexRange[0]);
}



FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
	WeakFormHMatAssemblyHelper);


}
