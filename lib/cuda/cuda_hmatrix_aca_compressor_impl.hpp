// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_CUDA_HMATRIX_ACA_COMPRESSOR_IMPL_HPP
#define HMAT_CUDA_HMATRIX_ACA_COMPRESSOR_IMPL_HPP

#include "cuda_hmatrix_aca_compressor.hpp"

#include "cuda_grid.hpp"

#include "../assembly/local_dof_lists_cache.hpp"

#include "../fiber/conjugate.hpp"

#include "../hmat/block_cluster_tree.hpp"
#include "../hmat/hmatrix_low_rank_data.hpp"

#include <tbb/task_group.h>
#include <tbb/concurrent_unordered_map.h>

#include <boost/numeric/conversion/cast.hpp>

#include <numeric>

// Helper functions and classes
namespace {

} // namespace

namespace hmat {

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::
CudaHMatrixAcaCompressor(
    const Bempp::Space<BasisFunctionType> &testSpace,
    const Bempp::Space<BasisFunctionType> &trialSpace,
    const shared_ptr<const Fiber::CollectionOfKernels<KernelType>> &kernel,
    const Fiber::Shapeset<BasisFunctionType> &testShapeset,
    const Fiber::Shapeset<BasisFunctionType> &trialShapeset,
    const Matrix<typename ScalarTraits<BasisFunctionType>::RealType>
    &localTestQuadPoints,
    const Matrix<typename ScalarTraits<BasisFunctionType>::RealType>
    &localTrialQuadPoints,
    const std::vector<typename ScalarTraits<BasisFunctionType>::RealType>
    &testQuadWeights,
    const std::vector<typename ScalarTraits<BasisFunctionType>::RealType>
    &trialQuadWeights,
    const shared_ptr<DefaultBlockClusterTreeType> blockClusterTree,
    const std::vector<LocalAssembler*> &assemblers,
    const std::vector<const DiscreteLinearOp*> &sparseTermsToAdd,
    const std::vector<ResultType> &denseTermsMultipliers,
    const std::vector<ResultType> &sparseTermsMultipliers,
    const double eps, const double zeroTol, const unsigned int maxRank,
    const DataAccessor<ResultType, N> &denseCompressor,
    const Bempp::CudaOptions &cudaOptions)
    : m_testSpace(testSpace), m_trialSpace(trialSpace),
      m_blockClusterTree(blockClusterTree), m_assemblers(assemblers),
      m_sparseTermsToAdd(sparseTermsToAdd),
      m_denseTermsMultipliers(denseTermsMultipliers),
      m_sparseTermsMultipliers(sparseTermsMultipliers),
      m_eps(eps), m_zeroTol(zeroTol), m_maxRank(maxRank),
      m_hMatrixDenseCompressor(denseCompressor),
      m_generalLocalTestDofCount(testShapeset.size()),
      m_testDofListsCache(new Bempp::LocalDofListsCache<BasisFunctionType>(
          m_testSpace,
          blockClusterTree->rowClusterTree()->hMatDofToOriginalDofMap(), true)),
      m_trialDofListsCache(new Bempp::LocalDofListsCache<BasisFunctionType>(
          m_trialSpace,
          blockClusterTree->columnClusterTree()->hMatDofToOriginalDofMap(),
          true)) {

  typedef typename ScalarTraits<CudaBasisFunctionType>::RealType CudaCoordinateType;

  for (size_t i = 0; i < assemblers.size(); ++i)
    if (!assemblers[i])
      throw std::invalid_argument(
          "CudaHMatrixAcaCompressor::CudaHMatrixAcaCompressor(): "
          "no elements of the 'assemblers' vector may be null");

  for (size_t i = 0; i < sparseTermsToAdd.size(); ++i)
    if (!sparseTermsToAdd[i])
      throw std::invalid_argument(
          "CudaHMatrixAcaCompressor::CudaHMatrixAcaCompressor(): "
          "no elements of the 'sparseTermsToAdd' vector may be null");

  m_accessedEntryCount = 0;

  const std::vector<int> &deviceIds = cudaOptions.devices();
  const unsigned int deviceCount = deviceIds.size();

  const unsigned int trialDofCount = trialShapeset.size();
  const unsigned int testDofCount = testShapeset.size();

  const unsigned int trialPointCount = localTrialQuadPoints.cols();
  const unsigned int testPointCount = localTestQuadPoints.cols();

  m_cudaIntegrators.resize(deviceCount);
  m_results.resize(deviceCount);

  // Get maximum number of element pairs which can be treated on the device
  // simultaneously
  m_maxActiveElemPairCount = cudaOptions.chunkElemPairCount();

  // Let the chunk size be a multiple of the warp size
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  const int warpSize = prop.warpSize;
  m_maxActiveElemPairCount =
      std::max((m_maxActiveElemPairCount / warpSize) * warpSize, size_t(warpSize));
  std::cout << "maxActiveElemPairCount = " << m_maxActiveElemPairCount << std::endl;

  // Loop over all devices
  for (int device = 0; device < deviceCount; ++device) {

    const int deviceId = deviceIds[device];
    cu_verify( cudaSetDevice(deviceId) );

    // Push raw grid data to the device
    shared_ptr<Bempp::CudaGrid<CudaCoordinateType>> trialGrid =
        trialSpace.grid()->template pushToDevice<CudaCoordinateType>(deviceId);
    shared_ptr<Bempp::CudaGrid<CudaCoordinateType>> testGrid =
        testSpace.grid()->template pushToDevice<CudaCoordinateType>(deviceId);

    // Allocate pinned host memory on the device
    size_t size = m_maxActiveElemPairCount
        * trialDofCount * testDofCount * sizeof(CudaResultType);
    if (Bempp::ScalarTraits<ResultType>::isComplex) size *= 2;
    cu_verify(
        cudaHostAlloc((void**)&(m_results[device]), size, cudaHostAllocMapped) );

    // Create CUDA integrator
    m_cudaIntegrators[device] = boost::make_shared<Fiber::CudaIntegrator<
        BasisFunctionType, KernelType, ResultType,
        CudaBasisFunctionType, CudaKernelType, CudaResultType>>(
            localTestQuadPoints, localTrialQuadPoints,
            testQuadWeights, trialQuadWeights,
            testShapeset, trialShapeset,
            testGrid, trialGrid,
            m_maxActiveElemPairCount,
            kernel,
            deviceId, cudaOptions);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
void CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::
compressAllBlocks(ParallelDataContainer &hMatrixData) {

  hMatrixData.clear();

  std::function<void(const node_t &node,
      std::vector<node_t> &admissibleLeafNodeList)> compressionInitFun =
      [&](const node_t &node,
          std::vector<node_t> &admissibleLeafNodeList) {

        if (node->isLeaf()) {

          shared_ptr<HMatrixData<ResultType>> nodeData;

          if (!(node->data().admissible)) {

            m_hMatrixDenseCompressor.compressBlock(*node, nodeData);

          } else {

            IndexRangeType rowClusterRange;
            IndexRangeType columnClusterRange;
            size_t numberOfRows;
            size_t numberOfColumns;

            getBlockClusterTreeNodeDimensions(*node, rowClusterRange,
                                              columnClusterRange, numberOfRows,
                                              numberOfColumns);

            nodeData.reset(new HMatrixLowRankData<ResultType>());
            Matrix<ResultType> &A =
                static_cast<HMatrixLowRankData<ResultType>*>(nodeData.get())->A();
            Matrix<ResultType> &B =
                static_cast<HMatrixLowRankData<ResultType>*>(nodeData.get())->B();

            A.resize(numberOfRows, 0);
            B.resize(0, numberOfColumns);

            // Add the current leaf to the list of admissible leafs
            admissibleLeafNodeList.push_back(node);

            m_rowApproxCounters.push_back(std::vector<size_t>(numberOfRows));
            m_colApproxCounters.push_back(std::vector<size_t>(numberOfColumns));
            m_maxIterations.push_back(std::min(static_cast<size_t>(m_maxRank),
                std::min(numberOfRows, numberOfColumns)));
          }

          hMatrixData[node] = nodeData;

        } else {

          // In case the current node is no leaf, loop over it's 4 childs
          for (int child = 0; child < 4; ++child) {
            compressionInitFun(node->child(child), admissibleLeafNodeList);
          }
        }
      };

  // Create a list of all admissible leaf nodes
  std::vector<node_t> admissibleLeafNodeList;
  admissibleLeafNodeList.clear();

  // Initialize the compression algorithm and compute the non-admissible node data
  compressionInitFun(m_blockClusterTree->root(), admissibleLeafNodeList);

  // Get the total number of adissible leaf nodes which have to be approximated
  const size_t admissibleLeafNodeCount = admissibleLeafNodeList.size();
  std::cout << "admissibleLeafNodeCount = " << admissibleLeafNodeCount << std::endl;

  // Create further vectors to keep track of the compression state of the
  // admissible leaf nodes
  m_origRows.resize(admissibleLeafNodeCount);
  m_origCols.resize(admissibleLeafNodeCount);
  m_blockNorms.resize(admissibleLeafNodeCount, 0.);
  m_finished.resize(admissibleLeafNodeCount, false);
  m_acaAlgorithmStates.resize(admissibleLeafNodeCount, AcaAlgorithmStateType::START);
  m_acaStatuses.resize(admissibleLeafNodeCount);
  m_nextPivots.resize(admissibleLeafNodeCount, 0);
  m_modes.resize(admissibleLeafNodeCount, ModeType::ROW);
  m_rowTrialPassed.resize(admissibleLeafNodeCount, false);
  m_columnTrialPassed.resize(admissibleLeafNodeCount, false);
  m_rowTrialCounts.resize(admissibleLeafNodeCount, 0);
  m_columnTrialCounts.resize(admissibleLeafNodeCount, 0);
  m_MAX_ROW_TRIAL_COUNTS.resize(admissibleLeafNodeCount, 1);
  m_MAX_COLUMN_TRIAL_COUNTS.resize(admissibleLeafNodeCount, 1);
  m_rowModeZeroTermination.resize(admissibleLeafNodeCount, false);
  m_columnModeZeroTermination.resize(admissibleLeafNodeCount, false);

  m_acaConverged.resize(admissibleLeafNodeCount, false);
  m_acaNextPivots.resize(admissibleLeafNodeCount, 0);
  m_acaIterationCounts.resize(admissibleLeafNodeCount, 0);
  m_acaCrossStatuses.resize(admissibleLeafNodeCount);
  m_rows.resize(admissibleLeafNodeCount);
  m_cols.resize(admissibleLeafNodeCount);

  // Create a list of leaf nodes whose element pairs are currently queued
  // for integral computation on the device
  std::vector<node_t> queuedLeafNodeList(admissibleLeafNodeCount);
  m_queuedNodalElemPairCounts.resize(admissibleLeafNodeCount);

  // Initialise the element pair queues
  m_queuedTestElemPairIndices.resize(m_maxActiveElemPairCount);
  m_queuedTrialElemPairIndices.resize(m_maxActiveElemPairCount);

  // Point to the next device
  m_nextDevice = 0;

  // Reset counters
  m_queuedElemPairCount = 0;
  m_queuedNodeCount = 0;

  // Run the algorithm until all admissible leaf nodes are finished
  size_t finishedAdmissibleNodeCount = 0;
//  while (finishedAdmissibleNodeCount != admissibleLeafNodeCount) {

    std::cout << "PART 1" << std::endl;

    // STEP 1: Go through all admissible leaf nodes in the list.
    for (size_t nd = 0; nd < admissibleLeafNodeCount; ++nd) {

      // In case this admissible leaf node is already finished, go to the next one
      if (m_finished[nd]) continue;

      const node_t &node = admissibleLeafNodeList[nd];

      if (m_maxIterations[nd] > 0) {

        const bool triggerIntegralComputation =
            addMatrixBlockIntegralsPart1(node, nd);

        if (!triggerIntegralComputation) {

          // Keep track of the nodes currently waiting for integral computation
          queuedLeafNodeList[m_queuedNodeCount] = node;
          m_queuedNodeCount++;

        } else {

          integrateAndAssembleQueuedMatrixBlocksPart1(
              admissibleLeafNodeList, queuedLeafNodeList);

          // Do this leaf node again, since it's integrals have not been considered
          nd--;
        }

      } else { // maxIterations <= 0

        setAndAnalyseAcaStatus(nd, finishedAdmissibleNodeCount);
      }
    }

    // Trigger integral computation one last time to empty the queue
    integrateAndAssembleQueuedMatrixBlocksPart1(
        admissibleLeafNodeList, queuedLeafNodeList);

    std::cout << "PART 2" << std::endl;

    // STEP 2: Go through all admissible leaf nodes in the list.
    for (size_t nd = 0; nd < admissibleLeafNodeCount; ++nd) {

      // In case this admissible leaf node is already finished, go to the next one
      if (m_finished[nd]) continue;

      const node_t &node = admissibleLeafNodeList[nd];

      const bool triggerIntegralComputation =
          addMatrixBlockIntegralsPart2(node, nd, hMatrixData);

      if (!triggerIntegralComputation) {

        // Find out whether the node is of ZERO cross status type.
        if (m_acaCrossStatuses[nd] == CrossStatusType::ZERO) {

          checkAcaConvergence(node, nd, hMatrixData, finishedAdmissibleNodeCount);

        } else {

          // Keep track of the nodes currently waiting for integral computation
          queuedLeafNodeList[m_queuedNodeCount] = node;
          m_queuedNodeCount++;
        }

      } else {

        integrateAndAssembleQueuedMatrixBlocksPart2(admissibleLeafNodeList,
            queuedLeafNodeList, hMatrixData, finishedAdmissibleNodeCount);

        // Do this leaf node again, since it's integrals have not been considered
        nd--;
      }
    }

    // Trigger integral computation one last time to empty the queue
    integrateAndAssembleQueuedMatrixBlocksPart2(admissibleLeafNodeList,
        queuedLeafNodeList, hMatrixData, finishedAdmissibleNodeCount);
//  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
bool CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::selectMinPivot(
    const Matrix<ResultType> &vec,
    const std::vector<std::size_t> &approximationCount, std::size_t &pivot,
    ModeType modus) const {

  int n = std::max(vec.rows(), vec.cols());
  pivot = n;
  int row = 0;
  int col = 0;
  double minVal = std::numeric_limits<double>::max();
  for (int i = 0; i < n; ++i) {
    if (modus == ModeType::ROW)
      col = i;
    else
      row = i;
    if (approximationCount[i] == 0 && std::abs(vec(row, col)) < minVal) {
      minVal = std::abs(vec(row, col));
      pivot = i;
    }
  }
  return (pivot != n);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
bool CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::checkConvergence(
    const Matrix<ResultType> &row, const Matrix<ResultType> &col,
    const double tol) const {

  return row.norm() * col.norm() < tol;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
double CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::
updateLowRankBlocksAndNorm(
    const Matrix<ResultType> &newRow, const Matrix<ResultType> &newCol,
    Matrix<ResultType> &A, Matrix<ResultType> &B,
    const double currentBlockNorm) const {

  // Compute norm update
  double newNormSquared = newRow.squaredNorm() * newCol.squaredNorm();

  if (A.cols() > 0 && B.rows() > 0) {
    double val1 =
        2 * (newRow * (B.adjoint() * (A.adjoint() * newCol))).real()(0, 0);
    double val2 = currentBlockNorm * currentBlockNorm;
    newNormSquared += val1 + val2;
  }

  // Now increase size of blocks
  Matrix<ResultType> Atmp(A.rows(), A.cols() + 1);
  Matrix<ResultType> Btmp(B.rows() + 1, B.cols());
  Atmp.leftCols(A.cols()) = A;
  Btmp.topRows(B.rows()) = B;
  A.swap(Atmp);
  B.swap(Btmp);

  A.rightCols(1) = newCol;
  B.bottomRows(1) = newRow;

  // Return 0 if due to rounding errors newNormSquared is smaller 0.
  return std::sqrt((newNormSquared > 0) ? newNormSquared : 0);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
bool CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::
addMatrixBlockIntegrals(
    const hmat::IndexRangeType &testIndexRange,
    const hmat::IndexRangeType &trialIndexRange,
    const hmat::DefaultBlockClusterTreeNodeType &blockClusterTreeNode) {

  bool triggerIntegralComputation = false;

  auto numberOfTestIndices = testIndexRange[1] - testIndexRange[0];
  auto numberOfTrialIndices = trialIndexRange[1] - trialIndexRange[0];

  shared_ptr<const Bempp::LocalDofLists<BasisFunctionType>> testDofLists =
      m_testDofListsCache->get(testIndexRange[0], numberOfTestIndices);
  shared_ptr<const Bempp::LocalDofLists<BasisFunctionType>> trialDofLists =
      m_trialDofListsCache->get(trialIndexRange[0], numberOfTrialIndices);

  // Necessary elements
  const std::vector<int> &testElementIndices = testDofLists->elementIndices;
  const std::vector<int> &trialElementIndices = trialDofLists->elementIndices;

//  if (blockClusterTreeNode.data().admissible) {
//    std::cout << "addMatrixBlockIntegrals" << std::endl;
//    for (int i = 0; i < testElementIndices.size(); ++i) {
//      std::cout << testElementIndices[i] << " " << std::flush;
//    }
//    std::cout << std::endl;
//    for (int i = 0; i < trialElementIndices.size(); ++i) {
//      std::cout << trialElementIndices[i] << " " << std::flush;
//    }
//    std::cout << std::endl;
//  }

  // Get the number of element pairs which would be added to the queue
  const size_t blockElemPairCount =
      testElementIndices.size() * trialElementIndices.size();

  // Check whether the maximum number of element pairs would be exceeded
  if ((m_queuedElemPairCount + blockElemPairCount) <= m_maxActiveElemPairCount) {

    m_queuedNodalElemPairCounts[m_queuedNodeCount] = blockElemPairCount;

    // Add element pairs to the queue
    for (size_t nTestElem = 0; nTestElem < testElementIndices.size();
         ++nTestElem) {
      const int activeTestElementIndex = testElementIndices[nTestElem];
      for (size_t nTrialElem = 0; nTrialElem < trialElementIndices.size();
           ++nTrialElem) {
        const int activeTrialElementIndex = trialElementIndices[nTrialElem];

        m_queuedTestElemPairIndices[m_queuedElemPairCount] = activeTestElementIndex;
        m_queuedTrialElemPairIndices[m_queuedElemPairCount] = activeTrialElementIndex;
        m_queuedElemPairCount++;
      }
    }

    triggerIntegralComputation = false;
    return triggerIntegralComputation;

  } else {

    triggerIntegralComputation = true;
    return triggerIntegralComputation;
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
void CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::
assembleMatrixBlock(
    const hmat::IndexRangeType &testIndexRange,
    const hmat::IndexRangeType &trialIndexRange,
    const hmat::DefaultBlockClusterTreeNodeType &blockClusterTreeNode,
    const size_t nd, Matrix<ResultType> &data) {

  auto numberOfTestIndices = testIndexRange[1] - testIndexRange[0];
  auto numberOfTrialIndices = trialIndexRange[1] - trialIndexRange[0];

  shared_ptr<const Bempp::LocalDofLists<BasisFunctionType>> testDofLists =
      m_testDofListsCache->get(testIndexRange[0], numberOfTestIndices);
  shared_ptr<const Bempp::LocalDofLists<BasisFunctionType>> trialDofLists =
      m_trialDofListsCache->get(trialIndexRange[0], numberOfTrialIndices);

  // Requested original matrix indices
  typedef typename Bempp::LocalDofLists<BasisFunctionType>::DofIndex DofIndex;
  const std::vector<DofIndex> &testOriginalIndices =
      testDofLists->originalIndices;
  const std::vector<DofIndex> &trialOriginalIndices =
      trialDofLists->originalIndices;

  // Necessary elements
  const std::vector<int> &testElementIndices = testDofLists->elementIndices;
  const std::vector<int> &trialElementIndices = trialDofLists->elementIndices;

//  if (blockClusterTreeNode.data().admissible) {
//    std::cout << "assembleMatrixBlock" << std::endl;
//    for (int i = 0; i < testElementIndices.size(); ++i) {
//      std::cout << testElementIndices[i] << " " << std::flush;
//    }
//    std::cout << std::endl;
//    for (int i = 0; i < trialElementIndices.size(); ++i) {
//      std::cout << trialElementIndices[i] << " " << std::flush;
//    }
//    std::cout << std::endl;
//  }

  // Necessary local dof indices in each element
  const std::vector<std::vector<Bempp::LocalDofIndex>> &testLocalDofs =
      testDofLists->localDofIndices;
  const std::vector<std::vector<Bempp::LocalDofIndex>> &trialLocalDofs =
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

  const size_t testElemCount = testElementIndices.size();
  const size_t trialElemCount = trialElementIndices.size();

  const size_t offset = std::accumulate(m_queuedNodalElemPairCounts.begin(),
      m_queuedNodalElemPairCounts.begin() + nd, 0);

  for (size_t testElem = 0; testElem < testElemCount; ++testElem) {
    const int activeTestElementIndex = testElementIndices[testElem];
    const size_t testDofCount = testLocalDofs[testElem].size();

    for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
      Bempp::LocalDofIndex activeTestLocalDof = testLocalDofs[testElem][testDof];
      BasisFunctionType activeTestLocalDofWeight =
          testLocalDofWeights[testElem][testDof];

      for (size_t trialElem = 0; trialElem < trialElemCount; ++trialElem) {
        const int activeTrialElementIndex = trialElementIndices[trialElem];
        const size_t trialDofCount = trialLocalDofs[trialElem].size();

        for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof) {

          const size_t index =
              trialLocalDofs[trialElem][trialDof] * m_generalLocalTestDofCount * m_queuedElemPairCount
              + testLocalDofs[testElem][testDof] * m_queuedElemPairCount
              + offset
              + trialElemCount * testElem
              + trialElem;

          data(blockRows[testElem][testDof],
               blockCols[trialElem][trialDof]) +=
//                m_denseTermsMultipliers[nTerm] *
              Fiber::conjugate(activeTestLocalDofWeight) *
              trialLocalDofWeights[trialElem][trialDof] *
              static_cast<ResultType>(m_results[m_nextDevice][index]);
        }
      }
    }
  }

//  if (blockClusterTreeNode.data().admissible) {
//    std::cout << "data = " << std::endl;
//    std::cout << data << std::endl;
//  }

  // Now, add the contributions of the sparse terms
  for (size_t nTerm = 0; nTerm < m_sparseTermsToAdd.size(); ++nTerm)
    m_sparseTermsToAdd[nTerm]->addBlock(
        // since m_indexWithGlobalDofs is set, these refer
        // to global DOFs
        testOriginalIndices, trialOriginalIndices,
        m_sparseTermsMultipliers[nTerm], data);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
void CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::
assembleComputedMatrixBlocksPart1(
    std::vector<node_t> &admissibleLeafNodeList,
    std::vector<node_t> &queuedLeafNodeList) {

  for (size_t nd = 0; nd < m_queuedNodeCount; ++nd) {

    const node_t &computedNode = queuedLeafNodeList[nd];

    std::ptrdiff_t index = std::find(admissibleLeafNodeList.begin(),
        admissibleLeafNodeList.end(), computedNode)
        - admissibleLeafNodeList.begin();

    IndexRangeType rowClusterRange;
    IndexRangeType columnClusterRange;
    size_t numberOfRows;
    size_t numberOfColumns;

    getBlockClusterTreeNodeDimensions(*computedNode, rowClusterRange,
                                      columnClusterRange, numberOfRows,
                                      numberOfColumns);
    size_t rowIndex;
    size_t columnIndex;

    IndexRangeType rowIndexRange;
    IndexRangeType columnIndexRange;

    if (m_modes[index] == ModeType::ROW) {

      rowIndex = rowClusterRange[0] + m_acaNextPivots[index];
      rowIndexRange = {{rowIndex, rowIndex + 1}};

      assembleMatrixBlock(rowIndexRange, columnClusterRange,
                          *computedNode, nd, m_origRows[index]);

    } else {

      columnIndex = columnClusterRange[0] + m_acaNextPivots[index];
      columnIndexRange = {{columnIndex, columnIndex + 1}};

      assembleMatrixBlock(rowClusterRange, columnIndexRange,
                          *computedNode, nd, m_origCols[index]);
    }
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
void CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::
assembleComputedMatrixBlocksPart2(
    std::vector<node_t> &admissibleLeafNodeList,
    std::vector<node_t> &queuedLeafNodeList,
    ParallelDataContainer &hMatrixData,
    size_t &finishedAdmissibleLeafNodeCount) {

  for (size_t nd = 0; nd < m_queuedNodeCount; ++nd) {

    const node_t &computedNode = queuedLeafNodeList[nd];

    std::ptrdiff_t index = std::find(admissibleLeafNodeList.begin(),
        admissibleLeafNodeList.end(), computedNode)
        - admissibleLeafNodeList.begin();

    IndexRangeType rowClusterRange;
    IndexRangeType columnClusterRange;
    size_t numberOfRows;
    size_t numberOfColumns;

    getBlockClusterTreeNodeDimensions(*computedNode, rowClusterRange,
                                      columnClusterRange, numberOfRows,
                                      numberOfColumns);

    Matrix<ResultType> &A = static_cast<HMatrixLowRankData<ResultType>*>(
        hMatrixData[computedNode].get())->A();
    Matrix<ResultType> &B = static_cast<HMatrixLowRankData<ResultType>*>(
        hMatrixData[computedNode].get())->B();

    size_t rowIndex;
    size_t columnIndex;

    std::ptrdiff_t maxRowInd;
    std::ptrdiff_t maxColInd;

    IndexRangeType rowIndexRange;
    IndexRangeType columnIndexRange;

    ResultType pivotValue;

    if (m_modes[index] == ModeType::ROW) {

      rowIndex = rowClusterRange[0] + m_acaNextPivots[index];

      m_rows[index].cwiseAbs().maxCoeff(&maxRowInd, &maxColInd);
      pivotValue = m_rows[index](0, maxColInd);
      columnIndex =
          boost::numeric_cast<size_t>(maxColInd) + columnClusterRange[0];
      columnIndexRange = {{columnIndex, columnIndex + 1}};

      assembleMatrixBlock(rowClusterRange, columnIndexRange,
                          *computedNode, nd, m_origCols[index]);

      m_cols[index] = m_origCols[index];
      if (A.cols() > 0 && B.rows() > 0)
        m_cols[index] -= A * B.col(columnIndex - columnClusterRange[0]);

      m_rows[index] /= pivotValue;

      // Find next pivot, make sure that current rowIndex
      // is associated with zero for the search.
      ResultType tmp = 0;
      std::swap(tmp, m_cols[index](rowIndex - rowClusterRange[0], 0));
      m_cols[index].cwiseAbs().maxCoeff(&maxRowInd, &maxColInd);
      std::swap(tmp, m_cols[index](rowIndex - rowClusterRange[0], 0));
      m_acaNextPivots[index] = boost::numeric_cast<size_t>(maxRowInd);

    } else {

      columnIndex = columnClusterRange[0] + m_acaNextPivots[index];

      m_cols[index].cwiseAbs().maxCoeff(&maxRowInd, &maxColInd);
      pivotValue = m_cols[index](maxRowInd, 0);
      rowIndex = boost::numeric_cast<size_t>(maxRowInd) + rowClusterRange[0];
      rowIndexRange = {{rowIndex, rowIndex + 1}};

      assembleMatrixBlock(rowIndexRange, columnClusterRange,
                          *computedNode, nd, m_origRows[index]);

      m_rows[index] = m_origRows[index];
      if (A.cols() > 0 && B.rows() > 0)
        m_rows[index] -= A.row(rowIndex - rowClusterRange[0]) * B;

      m_cols[index] /= pivotValue;

      // Find next pivot, make sure that current columnIndex
      // is associated with zero for the search.
      ResultType tmp = 0;
      std::swap(tmp, m_rows[index](0, columnIndex - columnClusterRange[0]));
      m_rows[index].cwiseAbs().maxCoeff(&maxRowInd, &maxColInd);
      std::swap(tmp, m_rows[index](0, columnIndex - columnClusterRange[0]));
      m_acaNextPivots[index] = boost::numeric_cast<size_t>(maxColInd);
    }

    m_rowApproxCounters[index][rowIndex - rowClusterRange[0]] += 1;
    m_colApproxCounters[index][columnIndex - columnClusterRange[0]] += 1;

    m_acaCrossStatuses[index] = CrossStatusType::SUCCESS;

    checkAcaConvergence(computedNode, index, hMatrixData, finishedAdmissibleLeafNodeCount);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
void CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::
integrateAndAssembleQueuedMatrixBlocksPart1(
    std::vector<node_t> &admissibleLeafNodeList,
    std::vector<node_t> &queuedLeafNodeList) {

  std::cout << "Trigger integral computation with "
            << m_queuedNodeCount << " queued nodes and "
            << m_queuedElemPairCount << " queued element pairs."
            << std::endl;

  // Send the element pair list to the device
  m_cudaIntegrators[m_nextDevice]->pushElemPairIndicesToDevice(
      m_queuedTestElemPairIndices, m_queuedTrialElemPairIndices,
      m_queuedElemPairCount);

  // Compute integrals on the device
  m_cudaIntegrators[m_nextDevice]->integrate(
      m_queuedElemPairCount, m_results[m_nextDevice]);

  // Assemble results
  assembleComputedMatrixBlocksPart1(
      admissibleLeafNodeList, queuedLeafNodeList);

  // Point to the next device
  m_nextDevice++;
  if (m_nextDevice >= m_cudaIntegrators.size()) m_nextDevice = 0;

  // Reset counters
  m_queuedElemPairCount = 0;
  m_queuedNodeCount = 0;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
void CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::
integrateAndAssembleQueuedMatrixBlocksPart2(
    std::vector<node_t> &admissibleLeafNodeList,
    std::vector<node_t> &queuedLeafNodeList,
    ParallelDataContainer &hMatrixData,
    size_t &finishedAdmissibleLeafNodeCount) {

  std::cout << "Trigger integral computation with "
            << m_queuedNodeCount << " queued nodes and "
            << m_queuedElemPairCount << " queued element pairs."
            << std::endl;

  // Send the element pair list to the device
  m_cudaIntegrators[m_nextDevice]->pushElemPairIndicesToDevice(
      m_queuedTestElemPairIndices, m_queuedTrialElemPairIndices,
      m_queuedElemPairCount);

  // Compute integrals on the device
  m_cudaIntegrators[m_nextDevice]->integrate(
      m_queuedElemPairCount, m_results[m_nextDevice]);

  // Assemble results
  assembleComputedMatrixBlocksPart2(admissibleLeafNodeList,
      queuedLeafNodeList, hMatrixData, finishedAdmissibleLeafNodeCount);

  // Point to the next device
  m_nextDevice++;
  if (m_nextDevice >= m_cudaIntegrators.size()) m_nextDevice = 0;

  // Reset counters
  m_queuedElemPairCount = 0;
  m_queuedNodeCount = 0;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
void CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::
setAndAnalyseAcaStatus(const size_t index,
                       size_t &finishedAdmissibleLeafNodeCount) {

  // Set the status of the ACA algorithm
  if (m_acaConverged[index]) {

    if (m_acaIterationCounts[index] == 0)
      m_acaStatuses[index] = AcaStatusType::CONVERGED_WITHOUT_ITERATION;
    else
      m_acaStatuses[index] = AcaStatusType::CONVERGED_AFTER_ITERATION;

  } else {
    m_acaStatuses[index] = AcaStatusType::RANK_LIMIT_REACHED;
  }

  // Analyse the status of the ACA algorithm and adjust states accordingly
  analyseAcaStatus(index, finishedAdmissibleLeafNodeCount);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
void CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::
checkAcaConvergence(const node_t &blockClusterTreeNode,
                    const size_t index,
                    ParallelDataContainer &hMatrixData,
                    size_t &finishedAdmissibleLeafNodeCount) {

  if (m_acaCrossStatuses[index] == CrossStatusType::ZERO) {

    if (m_acaIterationCounts[index] == 0)
     m_acaStatuses[index] = AcaStatusType::ZERO_TERMINATION_WITHOUT_ITERATION;
    else
     m_acaStatuses[index] = AcaStatusType::CONVERGED_AFTER_ITERATION;

    analyseAcaStatus(index, finishedAdmissibleLeafNodeCount);

  } else {

    m_acaConverged[index] =
        checkConvergence(m_rows[index], m_cols[index], m_eps * m_blockNorms[index]);

    if (m_acaConverged[index]) {

      setAndAnalyseAcaStatus(index, finishedAdmissibleLeafNodeCount);

    } else {

      Matrix<ResultType> &A = static_cast<HMatrixLowRankData<ResultType>*>(
          hMatrixData[blockClusterTreeNode].get())->A();
      Matrix<ResultType> &B = static_cast<HMatrixLowRankData<ResultType>*>(
          hMatrixData[blockClusterTreeNode].get())->B();

      m_blockNorms[index] = updateLowRankBlocksAndNorm(
          m_rows[index], m_cols[index], A, B, m_blockNorms[index]);

      m_maxIterations[index]--;
      m_acaIterationCounts[index]++;

      if (m_maxIterations[index] <= 0)
        setAndAnalyseAcaStatus(index, finishedAdmissibleLeafNodeCount);
    }
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
void CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::
analyseAcaStatus(const size_t index, size_t &finishedAdmissibleLeafNodeCount) {

  // Now test the status for different cases
  if (m_acaStatuses[index] == AcaStatusType::RANK_LIMIT_REACHED) {
    m_finished[index] = true; // Approximation does not work. So just stop.
    finishedAdmissibleLeafNodeCount++;
  }
  // The next should only happen if we are in ROW_TRIAL or COLUMN_TRIAL mode
  if (m_acaStatuses[index] == AcaStatusType::CONVERGED_WITHOUT_ITERATION) {
    if (m_acaAlgorithmStates[index] == AcaAlgorithmStateType::START)
      m_acaAlgorithmStates[index] = AcaAlgorithmStateType::ROW_TRIAL;
    else if (m_acaAlgorithmStates[index] == AcaAlgorithmStateType::ROW_TRIAL) {
      m_rowTrialPassed[index] = true;
      m_acaAlgorithmStates[index] = AcaAlgorithmStateType::COLUMN_TRIAL;
    } else if (m_acaAlgorithmStates[index] == AcaAlgorithmStateType::COLUMN_TRIAL) {
      m_columnTrialPassed[index] = true;
      m_acaAlgorithmStates[index] = AcaAlgorithmStateType::ROW_TRIAL;
    }
    // If both, the row trial and column trial are passed, finish.
    if (m_rowTrialPassed[index] && m_columnTrialPassed[index])
      m_finished[index] = true;
    finishedAdmissibleLeafNodeCount++;
  }
  // The next is the normal case after running ACA
  if (m_acaStatuses[index] == AcaStatusType::CONVERGED_AFTER_ITERATION ||
      m_acaStatuses[index] == AcaStatusType::ZERO_TERMINATION_AFTER_ITERATION) {
    // Have to rexecute row trial and column trial
    m_rowTrialPassed[index] = false;
    m_columnTrialPassed[index] = false;
    m_rowModeZeroTermination[index] = false;
    m_columnModeZeroTermination[index] = false;
    // Choose which trial to go into
    m_acaAlgorithmStates[index] = (m_modes[index] == ModeType::ROW)
        ? AcaAlgorithmStateType::ROW_TRIAL
        : AcaAlgorithmStateType::COLUMN_TRIAL;
  }
  if (m_acaStatuses[index] == AcaStatusType::ZERO_TERMINATION_WITHOUT_ITERATION) {
    if (m_modes[index] == ModeType::ROW) {
      m_rowModeZeroTermination[index] = true;
      m_acaAlgorithmStates[index] = AcaAlgorithmStateType::ROW_TRIAL;
    } else {
      m_columnModeZeroTermination[index] = true;
      m_acaAlgorithmStates[index] = AcaAlgorithmStateType::COLUMN_TRIAL;
    }
    if (m_rowModeZeroTermination[index] && m_columnModeZeroTermination[index])
      m_finished[index] = true;
      finishedAdmissibleLeafNodeCount++;
  }

  if (!m_finished[index] && m_acaAlgorithmStates[index] == AcaAlgorithmStateType::ROW_TRIAL) {
    if ((m_rowTrialCounts[index] < m_MAX_ROW_TRIAL_COUNTS[index]) &&
        selectMinPivot(m_origRows[index], m_colApproxCounters[index], m_nextPivots[index], ModeType::ROW)) {
      m_modes[index] = ModeType::COL;
      m_rowTrialCounts[index]++;
    } else {
      // All columns have been approximated
      m_finished[index] = true;
      finishedAdmissibleLeafNodeCount++;
    }
  }
  if (!m_finished[index] && m_acaAlgorithmStates[index] == AcaAlgorithmStateType::COLUMN_TRIAL) {
    if ((m_columnTrialCounts[index] < m_MAX_COLUMN_TRIAL_COUNTS[index]) &&
        selectMinPivot(m_origCols[index], m_rowApproxCounters[index], m_nextPivots[index], ModeType::COL)) {
      m_modes[index] = ModeType::ROW;
      m_columnTrialCounts[index]++;
    } else {
      // All columns have been approximated
      m_finished[index] = true;
      finishedAdmissibleLeafNodeCount++;
    }
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
bool CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::
addMatrixBlockIntegralsPart1(const node_t &blockClusterTreeNode,
                             const size_t index) {

  const auto &rowClusterRange =
      blockClusterTreeNode->data().rowClusterTreeNode->data().indexRange;

  const auto &columnClusterRange =
      blockClusterTreeNode->data().columnClusterTreeNode->data().indexRange;

  size_t rowIndex;
  size_t columnIndex;

  IndexRangeType rowIndexRange;
  IndexRangeType columnIndexRange;

  if (m_modes[index] == ModeType::ROW) {

    rowIndex = rowClusterRange[0] + m_acaNextPivots[index];
    rowIndexRange = {{rowIndex, rowIndex + 1}};

    return addMatrixBlockIntegrals(
        rowIndexRange, columnClusterRange, *blockClusterTreeNode);

  } else {

    columnIndex = columnClusterRange[0] + m_acaNextPivots[index];
    columnIndexRange = {{columnIndex, columnIndex + 1}};

    return addMatrixBlockIntegrals(
        rowClusterRange, columnIndexRange, *blockClusterTreeNode);
  }
}

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
bool CudaHMatrixAcaCompressor<BasisFunctionType, KernelType, ResultType,
CudaBasisFunctionType, CudaKernelType, CudaResultType, N>::
addMatrixBlockIntegralsPart2(const node_t &blockClusterTreeNode,
                             const size_t index,
                             ParallelDataContainer &hMatrixData) {

  Matrix<ResultType> &A =
      static_cast<HMatrixLowRankData<ResultType>*>(hMatrixData[blockClusterTreeNode].get())->A();
  Matrix<ResultType> &B =
      static_cast<HMatrixLowRankData<ResultType>*>(hMatrixData[blockClusterTreeNode].get())->B();

  const auto &rowClusterRange =
      blockClusterTreeNode->data().rowClusterTreeNode->data().indexRange;

  const auto &columnClusterRange =
      blockClusterTreeNode->data().columnClusterTreeNode->data().indexRange;

  size_t rowIndex;
  size_t columnIndex;

  std::ptrdiff_t maxRowInd;
  std::ptrdiff_t maxColInd;

  IndexRangeType rowIndexRange;
  IndexRangeType columnIndexRange;

  if (m_modes[index] == ModeType::ROW) {

    rowIndex = rowClusterRange[0] + m_acaNextPivots[index];

    m_rows[index] = m_origRows[index];
    if (A.cols() > 0 && B.rows() > 0)
      m_rows[index] -= A.row(rowIndex - rowClusterRange[0]) * B;

    if (m_rows[index].norm() <= m_zeroTol) {

      m_acaCrossStatuses[index] = CrossStatusType::ZERO;

      // Return at this point and don't add integrals
      return false;

    } else {
      m_rows[index].cwiseAbs().maxCoeff(&maxRowInd, &maxColInd);
      columnIndex =
          boost::numeric_cast<size_t>(maxColInd) + columnClusterRange[0];
      columnIndexRange = {{columnIndex, columnIndex + 1}};

      return addMatrixBlockIntegrals(
          rowClusterRange, columnIndexRange, *blockClusterTreeNode);
    }

  } else {

    columnIndex = columnClusterRange[0] + m_acaNextPivots[index];

    m_cols[index] = m_origCols[index];
    if (A.cols() > 0 && B.rows() > 0)
      m_cols[index] -= A * B.col(columnIndex - columnClusterRange[0]);

    if (m_cols[index].norm() <= m_zeroTol) {

      m_acaCrossStatuses[index] = CrossStatusType::ZERO;
      // Return at this point and don't add integrals
      return false;

    } else {

      m_cols[index].cwiseAbs().maxCoeff(&maxRowInd, &maxColInd);
      rowIndex = boost::numeric_cast<size_t>(maxRowInd) + rowClusterRange[0];
      rowIndexRange = {{rowIndex, rowIndex + 1}};

      return addMatrixBlockIntegrals(
          rowIndexRange, columnClusterRange, *blockClusterTreeNode);
    }
  }
}

}

#endif
