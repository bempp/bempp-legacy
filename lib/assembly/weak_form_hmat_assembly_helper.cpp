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
#include "local_dof_lists_cache.hpp"
#include "discrete_boundary_operator.hpp"

#include "../common/eigen_support.hpp"

#include "../hmat/block_cluster_tree.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/types.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/conjugate.hpp"

#include "../cuda/cuda_default_local_assembler_for_integral_operators_on_surfaces.hpp"

namespace Bempp {

using Fiber::conjugate;
}

namespace Bempp {

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
WeakFormHMatAssemblyHelper<BasisFunctionType, KernelType, ResultType,
    CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    WeakFormHMatAssemblyHelper(
        const Space<BasisFunctionType> &testSpace,
        const Space<BasisFunctionType> &trialSpace,
        const shared_ptr<const Fiber::CollectionOfKernels<KernelType>> &kernel,
        const Fiber::Shapeset<BasisFunctionType> &testShapeset,
        const Fiber::Shapeset<BasisFunctionType> &trialShapeset,
        const shared_ptr<hmat::DefaultBlockClusterTreeType> blockClusterTree,
        const std::vector<LocalAssembler *> &assemblers,
        const std::vector<const DiscreteLinOp *> &sparseTermsToAdd,
        const std::vector<ResultType> &denseTermsMultipliers,
        const std::vector<ResultType> &sparseTermsMultipliers,
        const Context<BasisFunctionType, ResultType> &context)
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
          true)),
      m_isCudaEnabled(context.assemblyOptions().isCudaEnabled()),
      m_cudaMinBlockSize(context.globalParameterList().template get<int>("options.hmat.cudaMinBlockSize")),
      m_deviceCount(context.cudaOptions().devices().size()), m_nextDevice(0),
      m_cpuBlockCount(0), m_accessedCpuEntryCount(0)/*,
      m_allocationTimer(0), m_integrationTimer(0), m_assemblyTimer(0)*/ {

  m_cudaBlockCount.resize(m_deviceCount, 0);
  m_accessedCudaEntryCount.resize(m_deviceCount, 0);

  for (size_t i = 0; i < assemblers.size(); ++i)
    if (!assemblers[i])
      throw std::invalid_argument(
          "WeakFormAcaAssemblyHelper::WeakFormHMatAssemblyHelper(): "
          "no elements of the 'assemblers' vector may be null");
  for (size_t i = 0; i < sparseTermsToAdd.size(); ++i)
    if (!sparseTermsToAdd[i])
      throw std::invalid_argument(
          "WeakFormAcaAssemblyHelper::WeakFormHMatAssemblyHelper(): "
          "no elements of the 'sparseTermsToAdd' vector may be null");
  m_accessedEntryCount = 0;

  if (m_isCudaEnabled) {
      m_cudaAssembler = new CudaLocalAssembler(
          testSpace, trialSpace, kernel, testShapeset, trialShapeset,
          context);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
WeakFormHMatAssemblyHelper<BasisFunctionType, KernelType, ResultType,
    CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    ~WeakFormHMatAssemblyHelper() {

  if (m_isCudaEnabled) {
    delete m_cudaAssembler;
    m_cudaAssembler = NULL;
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
typename WeakFormHMatAssemblyHelper<BasisFunctionType, KernelType, ResultType,
    CudaBasisFunctionType, CudaKernelType, CudaResultType>::MagnitudeType
WeakFormHMatAssemblyHelper<BasisFunctionType, KernelType, ResultType,
    CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    estimateMinimumDistance(const hmat::DefaultBlockClusterTreeNodeType
                                &blockClusterTreeNode) const {

  MagnitudeType dist =
      MagnitudeType(blockClusterTreeNode.data()
                        .rowClusterTreeNode->data()
                        .boundingBox.distance(blockClusterTreeNode.data()
                                                  .columnClusterTreeNode->data()
                                                  .boundingBox));
  return dist;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
double WeakFormHMatAssemblyHelper<BasisFunctionType, KernelType, ResultType,
    CudaBasisFunctionType, CudaKernelType, CudaResultType>::scale(
    const hmat::DefaultBlockClusterTreeNodeType &blockClusterTreeNode) const {

  MagnitudeType result = 0;
  MagnitudeType dist = this->estimateMinimumDistance(blockClusterTreeNode);

  for (size_t nTerm = 0; nTerm < m_assemblers.size(); ++nTerm)
    result = std::max(result, m_assemblers[nTerm]->estimateRelativeScale(dist));
  return static_cast<double>(result);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
void WeakFormHMatAssemblyHelper<BasisFunctionType, KernelType, ResultType,
    CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    computeMatrixBlock(
        const hmat::IndexRangeType &testIndexRange,
        const hmat::IndexRangeType &trialIndexRange,
        const hmat::DefaultBlockClusterTreeNodeType &blockClusterTreeNode,
        Matrix<ResultType> &data) const {

  auto numberOfTestIndices = testIndexRange[1] - testIndexRange[0];
  auto numberOfTrialIndices = trialIndexRange[1] - trialIndexRange[0];
//  printf("numberOfTestIndices = %d\n" , numberOfTestIndices);
//  printf("numberOfTrialIndices = %d\n", numberOfTrialIndices);

  m_accessedEntryCount += numberOfTestIndices * numberOfTrialIndices;

  data.resize(numberOfTestIndices, numberOfTrialIndices);
  data.setZero();

  const size_t cudaBlockSize =
//      testElementIndices.size() * trialElementIndices.size();
      numberOfTestIndices * numberOfTrialIndices;

  if (blockClusterTreeNode.data().admissible && m_isCudaEnabled &&
      cudaBlockSize >= m_cudaMinBlockSize) { // Use CUDA

    // Requested original matrix indices
    typedef typename LocalDofLists<BasisFunctionType>::DofIndex DofIndex;
    const std::vector<DofIndex> &testOriginalIndices =
        m_testDofListsCache->getOriginalIndices(testIndexRange[0], numberOfTestIndices);
    const std::vector<DofIndex> &trialOriginalIndices =
        m_trialDofListsCache->getOriginalIndices(trialIndexRange[0], numberOfTrialIndices);
  //    printf(" testOriginalIndices(%d) = ",  testOriginalIndices.size()); for (int i = 0; i <  testOriginalIndices.size(); ++i) { printf("%d ",  testOriginalIndices[i]); } printf("\n"); fflush(stdout);
  //    printf("trialOriginalIndices(%d) = ", trialOriginalIndices.size()); for (int i = 0; i < trialOriginalIndices.size(); ++i) { printf("%d ", trialOriginalIndices[i]); } printf("\n"); fflush(stdout);

    std::vector<std::vector<LocalDof>> testRawLocalDofs;
    std::vector<std::vector<BasisFunctionType>> testRawLocalDofWeights;
    m_testSpace.global2localDofs(testOriginalIndices, testRawLocalDofs, testRawLocalDofWeights);
//    printf("testRawLocalDofs(%d) = \n", testRawLocalDofs.size()); for (int i = 0; i < testRawLocalDofs.size(); ++i) { for (int j = 0; j < testRawLocalDofs[i].size(); ++j) { printf("[%d, %d] ", testRawLocalDofs[i][j].entityIndex, testRawLocalDofs[i][j].dofIndex); } printf("\n"); } printf("\n"); fflush(stdout);

    std::vector<std::vector<LocalDof>> trialRawLocalDofs;
    std::vector<std::vector<BasisFunctionType>> trialRawLocalDofWeights;
    m_trialSpace.global2localDofs(trialOriginalIndices, trialRawLocalDofs, trialRawLocalDofWeights);
//    printf("trialRawLocalDofs(%d) = \n", trialRawLocalDofs.size()); for (int i = 0; i < trialRawLocalDofs.size(); ++i) { for (int j = 0; j < trialRawLocalDofs[i].size(); ++j) { printf("[%d, %d] ", trialRawLocalDofs[i][j].entityIndex, trialRawLocalDofs[i][j].dofIndex); } printf("\n"); } printf("\n"); fflush(stdout);

    unsigned int nextDevice;
    {
      MutexType::scoped_lock lock(m_deviceMutex);
      nextDevice = m_nextDevice;
      m_nextDevice++;
      if (m_nextDevice >= m_deviceCount) m_nextDevice = 0;
    }
//    m_cudaAssembler->evaluateLocalWeakForms(
//        testElementIndices, trialElementIndices, testLocalDofs, trialLocalDofs,
//        testLocalDofWeights, trialLocalDofWeights, blockRows, blockCols,
//        data, nextDevice/*, m_allocationTimer, m_integrationTimer, m_assemblyTimer*/);
    m_cudaAssembler->evaluateLocalWeakForms(testRawLocalDofs, trialRawLocalDofs, data, nextDevice);

    m_cudaBlockCount[nextDevice]++;
    m_accessedCudaEntryCount[nextDevice] += numberOfTestIndices * numberOfTrialIndices;

  } else { // Do not use CUDA

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
  //    printf(" testOriginalIndices(%d) = ",  testOriginalIndices.size()); for (int i = 0; i <  testOriginalIndices.size(); ++i) { printf("%d ",  testOriginalIndices[i]); } printf("\n"); fflush(stdout);
  //    printf("trialOriginalIndices(%d) = ", trialOriginalIndices.size()); for (int i = 0; i < trialOriginalIndices.size(); ++i) { printf("%d ", trialOriginalIndices[i]); } printf("\n"); fflush(stdout);

    // Necessary elements
    const std::vector<int> & testElementIndices =  testDofLists->elementIndices;
    const std::vector<int> &trialElementIndices = trialDofLists->elementIndices;
//    printf(" testElemIndices(%d) = ",  testElementIndices.size()); for (int i = 0; i <  testElementIndices.size(); ++i) { printf("%d ",  testElementIndices[i]); } printf("\n"); fflush(stdout);
//    printf("trialElemIndices(%d) = ", trialElementIndices.size()); for (int i = 0; i < trialElementIndices.size(); ++i) { printf("%d ", trialElementIndices[i]); } printf("\n"); fflush(stdout);

    // Necessary local dof indices in each element
    const std::vector<std::vector<LocalDofIndex>> &testLocalDofs =
        testDofLists->localDofIndices;
    const std::vector<std::vector<LocalDofIndex>> &trialLocalDofs =
        trialDofLists->localDofIndices;
//    printf("testLocalDofs(%d) = \n",   testLocalDofs.size()); for (int i = 0; i <  testLocalDofs.size(); ++i) { for (int j = 0; j <  testLocalDofs[i].size(); ++j) { printf("%d ",  testLocalDofs[i][j]); } printf("\n"); } printf("\n"); fflush(stdout);
//    printf("trialLocalDofs(%d) = \n", trialLocalDofs.size()); for (int i = 0; i < trialLocalDofs.size(); ++i) { for (int j = 0; j < trialLocalDofs[i].size(); ++j) { printf("%d ", trialLocalDofs[i][j]); } printf("\n"); } printf("\n"); fflush(stdout);

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
    const std::vector<std::vector<int>> &blockRows =  testDofLists->arrayIndices;
    const std::vector<std::vector<int>> &blockCols = trialDofLists->arrayIndices;
//    printf("blockRows(%d) = \n", blockRows.size()); for (int i = 0; i < blockRows.size(); ++i) { for (int j = 0; j < blockRows[i].size(); ++j) { printf("%d ", blockRows[i][j]); } printf("\n"); } printf("\n"); fflush(stdout);
//    printf("blockCols(%d) = \n", blockCols.size()); for (int i = 0; i < blockCols.size(); ++i) { for (int j = 0; j < blockCols[i].size(); ++j) { printf("%d ", blockCols[i][j]); } printf("\n"); } printf("\n"); fflush(stdout);

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

    m_cpuBlockCount++;
    m_accessedCpuEntryCount += numberOfTestIndices * numberOfTrialIndices;

    // Now, add the contributions of the sparse terms
    for (size_t nTerm = 0; nTerm < m_sparseTermsToAdd.size(); ++nTerm)
      m_sparseTermsToAdd[nTerm]->addBlock(
          // since m_indexWithGlobalDofs is set, these refer
          // to global DOFs
          testOriginalIndices, trialOriginalIndices,
          m_sparseTermsMultipliers[nTerm], data);
  }

}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType>
void WeakFormHMatAssemblyHelper<BasisFunctionType, KernelType, ResultType,
    CudaBasisFunctionType, CudaKernelType, CudaResultType>::
    getStatistics(
        std::vector<size_t> &cudaBlockCount, size_t &cpuBlockCount,
        std::vector<size_t> &accessedCudaEntryCount, size_t &accessedCpuEntryCount,
        double &allocationTimer, double &integrationTimer, double &assemblyTimer) const {

    size_t totalCudaBlockCount = 0;
    cudaBlockCount.resize(m_deviceCount); accessedCudaEntryCount.resize(m_deviceCount);
    for (int device = 0; device < m_deviceCount; ++device) {
              totalCudaBlockCount   +=         m_cudaBlockCount[device];
              cudaBlockCount[device] =         m_cudaBlockCount[device];
      accessedCudaEntryCount[device] = m_accessedCudaEntryCount[device];
    }
           cpuBlockCount =         m_cpuBlockCount;
   accessedCpuEntryCount = m_accessedCpuEntryCount;

//    allocationTimer =  std::chrono::duration_cast<std::chrono::seconds>( m_allocationTimer).count() / totalCudaBlockCount;
//   integrationTimer =  std::chrono::duration_cast<std::chrono::seconds>(m_integrationTimer).count() / totalCudaBlockCount;
//      assemblyTimer =  std::chrono::duration_cast<std::chrono::seconds>(   m_assemblyTimer).count() / totalCudaBlockCount;
}

// Explicit instantiations
template class WeakFormHMatAssemblyHelper<float, float, float, float, float, float>;
template class WeakFormHMatAssemblyHelper<float, float, std::complex<float>, float, float, float>;
template class WeakFormHMatAssemblyHelper<float, std::complex<float>, std::complex<float>, float, float, float>;
template class WeakFormHMatAssemblyHelper<std::complex<float>, float, std::complex<float>, float, float, float>;
template class WeakFormHMatAssemblyHelper<std::complex<float>, std::complex<float>, std::complex<float>, float, float, float>;

template class WeakFormHMatAssemblyHelper<double, double, double, double, double, double>;
template class WeakFormHMatAssemblyHelper<double, double, std::complex<double>, double, double, double>;
template class WeakFormHMatAssemblyHelper<double, std::complex<double>, std::complex<double>, double, double, double>;
template class WeakFormHMatAssemblyHelper<std::complex<double>, double, std::complex<double>, double, double, double>;
template class WeakFormHMatAssemblyHelper<std::complex<double>, std::complex<double>, std::complex<double>, double, double, double>;

template class WeakFormHMatAssemblyHelper<float, float, float, double, double, double>;
template class WeakFormHMatAssemblyHelper<float, float, std::complex<float>, double, double, double>;
template class WeakFormHMatAssemblyHelper<float, std::complex<float>, std::complex<float>, double, double, double>;
template class WeakFormHMatAssemblyHelper<std::complex<float>, float, std::complex<float>, double, double, double>;
template class WeakFormHMatAssemblyHelper<std::complex<float>, std::complex<float>, std::complex<float>, double, double, double>;

template class WeakFormHMatAssemblyHelper<double, double, double, float, float, float>;
template class WeakFormHMatAssemblyHelper<double, double, std::complex<double>, float, float, float>;
template class WeakFormHMatAssemblyHelper<double, std::complex<double>, std::complex<double>, float, float, float>;
template class WeakFormHMatAssemblyHelper<std::complex<double>, double, std::complex<double>, float, float, float>;
template class WeakFormHMatAssemblyHelper<std::complex<double>, std::complex<double>, std::complex<double>, float, float, float>;

}
