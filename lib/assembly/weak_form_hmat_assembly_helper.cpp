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
#include "discrete_boundary_operator.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/conjugate.hpp"
#include "../common/eigen_support.hpp"

namespace Bempp {

using Fiber::conjugate;
}

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
        const std::vector<ResultType> &sparseTermsMultipliers)
    : m_testSpace(testSpace), m_trialSpace(trialSpace),
      m_blockClusterTree(blockClusterTree), m_assemblers(assemblers),
      m_sparseTermsToAdd(sparseTermsToAdd),
      m_denseTermsMultipliers(denseTermsMultipliers),
      m_sparseTermsMultipliers(sparseTermsMultipliers),
      m_testDofListsCache(new LocalDofListsCache<BasisFunctionType>(
          m_testSpace,
          blockClusterTree->rowClusterTree()->hMatDofToOriginalDofMap(), true)),
      m_trialDofListsCache(new LocalDofListsCache<BasisFunctionType>(
          m_trialSpace,
          blockClusterTree->columnClusterTree()->hMatDofToOriginalDofMap(),
          true)) {

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
              .rowClusterTreeNode->data().boundingBox.distance(
              blockClusterTreeNode.data().
              columnClusterTreeNode->data().boundingBox));
}

template <typename BasisFunctionType, typename ResultType>
void
WeakFormHMatAssemblyHelper<BasisFunctionType, ResultType>::computeMatrixBlock(
    const hmat::IndexRangeType &testIndexRange,
    const hmat::IndexRangeType &trialIndexRange,
    const hmat::DefaultBlockClusterTreeNodeType &blockClusterTreeNode,
    Matrix<ResultType> &data) const {

  auto numberOfTestIndices = testIndexRange[1] - testIndexRange[0];
  auto numberOfTrialIndices = trialIndexRange[1] - trialIndexRange[0];

  m_accessedEntryCount += numberOfTestIndices * numberOfTrialIndices;

  const CoordinateType minDist = estimateMinimumDistance(blockClusterTreeNode);

  shared_ptr<const LocalDofLists<BasisFunctionType>> testDofLists =
      m_testDofListsCache->get(testIndexRange[0], numberOfTestIndices);
  shared_ptr<const LocalDofLists<BasisFunctionType>> trialDofLists =
      m_trialDofListsCache->get(trialIndexRange[0], numberOfTrialIndices);

  // Requested original matrix indices
  typedef typename LocalDofLists<BasisFunctionType>::DofIndex DofIndex;
  const std::vector<DofIndex> &testOriginalIndices =
      testDofLists->originalIndices;
  const std::vector<DofIndex> &trialOriginalIndices =
      trialDofLists->originalIndices;
  // Necessary elements
  const std::vector<int> &testElementIndices = testDofLists->elementIndices;
  const std::vector<int> &trialElementIndices = trialDofLists->elementIndices;
  // Necessary local dof indices in each element
  const std::vector<std::vector<LocalDofIndex>> &testLocalDofs =
      testDofLists->localDofIndices;
  const std::vector<std::vector<LocalDofIndex>> &trialLocalDofs =
      trialDofLists->localDofIndices;
  // Weights of local dofs in each element
  const std::vector<std::vector<BasisFunctionType>> &testLocalDofWeights =
      testDofLists->localDofWeights;
  const std::vector<std::vector<BasisFunctionType>> &trialLocalDofWeights =
      trialDofLists->localDofWeights;
  for (size_t i = 0; i < testLocalDofWeights.size(); ++i)
    for (size_t j = 0; j < testLocalDofWeights[i].size(); ++j)
      assert(std::abs(testLocalDofWeights[i][j]) > 0.);
  for (size_t i = 0; i < trialLocalDofWeights.size(); ++i)
    for (size_t j = 0; j < trialLocalDofWeights[i].size(); ++j)
      assert(std::abs(trialLocalDofWeights[i][j]) > 0.);

  // Corresponding row and column indices in the matrix to be calculated
  const std::vector<std::vector<int>> &blockRows = testDofLists->arrayIndices;
  const std::vector<std::vector<int>> &blockCols = trialDofLists->arrayIndices;

  data.resize(numberOfTestIndices, numberOfTrialIndices);
  data.setZero();

  // First, evaluate the contributions of the dense terms
  if (numberOfTrialIndices == 1) {
    // Only one column of the block needed. This means that we need only
    // one local DOF from just one or a few trialElements. Evaluate the
    // local weak form for one local trial DOF at a time.

    std::vector<Matrix<ResultType>> localResult;
    for (size_t nTrialElem = 0; nTrialElem < trialElementIndices.size();
         ++nTrialElem) {
      const int activeTrialElementIndex = trialElementIndices[nTrialElem];

      // The body of this loop will very probably only run once (single
      // local DOF per trial element)
      for (size_t nTrialDof = 0; nTrialDof < trialLocalDofs[nTrialElem].size();
           ++nTrialDof) {
        LocalDofIndex activeTrialLocalDof =
            trialLocalDofs[nTrialElem][nTrialDof];
        BasisFunctionType activeTrialLocalDofWeight =
            trialLocalDofWeights[nTrialElem][nTrialDof];
        for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm) {
          m_assemblers[nTerm]->evaluateLocalWeakForms(
              Fiber::TEST_TRIAL, testElementIndices, activeTrialElementIndex,
              activeTrialLocalDof, localResult, minDist);
          for (size_t nTestElem = 0; nTestElem < testElementIndices.size();
               ++nTestElem)
            for (size_t nTestDof = 0;
                 nTestDof < testLocalDofs[nTestElem].size(); ++nTestDof)
              data(blockRows[nTestElem][nTestDof], 0) +=
                  m_denseTermsMultipliers[nTerm] *
                  conjugate(testLocalDofWeights[nTestElem][nTestDof]) *
                  activeTrialLocalDofWeight *
                  localResult[nTestElem](testLocalDofs[nTestElem][nTestDof]);
        }
      }
    }
  } else if (numberOfTestIndices == 1) // very few testElements
  {
    // Only one row of the block needed. This means that we need only
    // one local DOF from just one or a few testElements. Evaluate the
    // local weak form for one local test DOF at a time.

    std::vector<Matrix<ResultType>> localResult;
    for (size_t nTestElem = 0; nTestElem < testElementIndices.size();
         ++nTestElem) {
      const int activeTestElementIndex = testElementIndices[nTestElem];
      // The body of this loop will very probably only run once (single
      // local DOF per test element)
      for (size_t nTestDof = 0; nTestDof < testLocalDofs[nTestElem].size();
           ++nTestDof) {
        LocalDofIndex activeTestLocalDof = testLocalDofs[nTestElem][nTestDof];
        BasisFunctionType activeTestLocalDofWeight =
            testLocalDofWeights[nTestElem][nTestDof];
        for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm) {
          m_assemblers[nTerm]->evaluateLocalWeakForms(
              Fiber::TRIAL_TEST, trialElementIndices, activeTestElementIndex,
              activeTestLocalDof, localResult, minDist);
          for (size_t nTrialElem = 0; nTrialElem < trialElementIndices.size();
               ++nTrialElem)
            for (size_t nTrialDof = 0;
                 nTrialDof < trialLocalDofs[nTrialElem].size(); ++nTrialDof)
              data(0, blockCols[nTrialElem][nTrialDof]) +=
                  m_denseTermsMultipliers[nTerm] *
                  conjugate(activeTestLocalDofWeight) *
                  trialLocalDofWeights[nTrialElem][nTrialDof] *
                  localResult[nTrialElem](
                      trialLocalDofs[nTrialElem][nTrialDof]);
        }
      }
    }
  } else if (numberOfTestIndices <= 32 &&
             numberOfTrialIndices <= 32) // a "fat" block
  {
    // The whole block or its submatrix needed. This means that we are
    // likely to need all or almost all local DOFs from most elements.
    // Evaluate the full local weak form for each pair of test and trial
    // elements and then select the entries that we need.

    Fiber::_2dArray<Matrix<ResultType>> localResult;
    for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm) {
      m_assemblers[nTerm]->evaluateLocalWeakForms(
          testElementIndices, trialElementIndices, localResult, minDist);
      for (size_t nTrialElem = 0; nTrialElem < trialElementIndices.size();
           ++nTrialElem)
        for (size_t nTrialDof = 0;
             nTrialDof < trialLocalDofs[nTrialElem].size(); ++nTrialDof)
          for (size_t nTestElem = 0; nTestElem < testElementIndices.size();
               ++nTestElem)
            for (size_t nTestDof = 0;
                 nTestDof < testLocalDofs[nTestElem].size(); ++nTestDof)
              data(blockRows[nTestElem][nTestDof],
                   blockCols[nTrialElem][nTrialDof]) +=
                  m_denseTermsMultipliers[nTerm] *
                  conjugate(testLocalDofWeights[nTestElem][nTestDof]) *
                  trialLocalDofWeights[nTrialElem][nTrialDof] *
                  localResult(nTestElem, nTrialElem)(
                      testLocalDofs[nTestElem][nTestDof],
                      trialLocalDofs[nTrialElem][nTrialDof]);
    }
  } else {
    std::vector<Matrix<ResultType>> localResult;
    for (size_t nTestElem = 0; nTestElem < testElementIndices.size();
         ++nTestElem) {
      const int activeTestElementIndex = testElementIndices[nTestElem];
      // The body of this loop will very probably only run once (single
      // local DOF per test element)
      for (size_t nTestDof = 0; nTestDof < testLocalDofs[nTestElem].size();
           ++nTestDof) {
        LocalDofIndex activeTestLocalDof = testLocalDofs[nTestElem][nTestDof];
        BasisFunctionType activeTestLocalDofWeight =
            testLocalDofWeights[nTestElem][nTestDof];
        for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm) {
          m_assemblers[nTerm]->evaluateLocalWeakForms(
              Fiber::TRIAL_TEST, trialElementIndices, activeTestElementIndex,
              activeTestLocalDof, localResult, minDist);
          for (size_t nTrialElem = 0; nTrialElem < trialElementIndices.size();
               ++nTrialElem)
            for (size_t nTrialDof = 0;
                 nTrialDof < trialLocalDofs[nTrialElem].size(); ++nTrialDof)
              data(blockRows[nTestElem][nTestDof],
                   blockCols[nTrialElem][nTrialDof]) +=
                  m_denseTermsMultipliers[nTerm] *
                  conjugate(activeTestLocalDofWeight) *
                  trialLocalDofWeights[nTrialElem][nTrialDof] *
                  localResult[nTrialElem](
                      trialLocalDofs[nTrialElem][nTrialDof]);
        }
      }
    }
  }

  // Now, add the contributions of the sparse terms
  for (size_t nTerm = 0; nTerm < m_sparseTermsToAdd.size(); ++nTerm)
    m_sparseTermsToAdd[nTerm]->addBlock(
        // since m_indexWithGlobalDofs is set, these refer
        // to global DOFs
        testOriginalIndices, trialOriginalIndices,
        m_sparseTermsMultipliers[nTerm], data);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    WeakFormHMatAssemblyHelper);
}
