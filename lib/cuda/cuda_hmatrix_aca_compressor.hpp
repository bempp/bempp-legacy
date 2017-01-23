// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_CUDA_HMATRIX_ACA_COMPRESSOR_HPP
#define HMAT_CUDA_HMATRIX_ACA_COMPRESSOR_HPP

#include "cuda_integrator.hpp"

#include "cuda_options.hpp"

#include "../hmat/common.hpp"
#include "../hmat/eigen_fwd.hpp"
#include "../hmat/block_cluster_tree.hpp"
#include "../hmat/hmatrix_data.hpp"
#include "../hmat/data_accessor.hpp"
#include "../hmat/hmatrix_dense_compressor.hpp"

#include <tbb/spin_mutex.h>

namespace Bempp {

/** \cond FORWARD_DECL */
template <typename ResultType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType> class LocalDofListsCache;
template <typename BasisFunctionType> class Space;
/** \endcond */

} // namespace Bempp

namespace hmat {

template <typename BasisFunctionType, typename KernelType, typename ResultType,
typename CudaBasisFunctionType, typename CudaKernelType, typename CudaResultType,
int N>
class CudaHMatrixAcaCompressor {
public:

  typedef Bempp::DiscreteBoundaryOperator<ResultType> DiscreteLinearOp;
  typedef Fiber::LocalAssemblerForIntegralOperators<ResultType> LocalAssembler;
  typedef shared_ptr<hmat::DefaultBlockClusterTreeNodeType> node_t;

  typedef tbb::concurrent_unordered_map<
      shared_ptr<BlockClusterTreeNode<N>>,
      shared_ptr<HMatrixData<ResultType>>,
      shared_ptr_hash<BlockClusterTreeNode<N>>>
          ParallelDataContainer;

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
      const Bempp::CudaOptions &cudaOptions);

  void compressAllBlocks(ParallelDataContainer &hMatrixData);

  enum class ModeType { ROW, COL };
  enum class CrossStatusType { SUCCESS, ZERO };
  enum class AcaStatusType {
    CONVERGED_AFTER_ITERATION,
    CONVERGED_WITHOUT_ITERATION,
    ZERO_TERMINATION_WITHOUT_ITERATION,
    ZERO_TERMINATION_AFTER_ITERATION,
    RANK_LIMIT_REACHED
  };
  enum class AcaAlgorithmStateType { START, ROW_TRIAL, COLUMN_TRIAL };

private:
  /** \cond PRIVATE */

  bool addMatrixBlockIntegrals(
      const hmat::IndexRangeType &testIndexRange,
      const hmat::IndexRangeType &trialIndexRange,
      const hmat::DefaultBlockClusterTreeNodeType &blockClusterTreeNode);

  void assembleMatrixBlock(
      const hmat::IndexRangeType &testIndexRange,
      const hmat::IndexRangeType &trialIndexRange,
      const hmat::DefaultBlockClusterTreeNodeType &blockClusterTreeNode,
      const size_t nd, Matrix<ResultType> &data);

  void assembleComputedMatrixBlocksPart1(
      std::vector<node_t> &admissibleLeafNodeList,
      std::vector<node_t> &queuedLeafNodeList);

  void assembleComputedMatrixBlocksPart2(
      std::vector<node_t> &admissibleLeafNodeList,
      std::vector<node_t> &queuedLeafNodeList,
      ParallelDataContainer &hMatrixData,
      size_t &finishedAdmissibleLeafNodeCount);

  void integrateAndAssembleQueuedMatrixBlocksPart1(
      std::vector<node_t> &admissibleLeafNodeList,
      std::vector<node_t> &queuedLeafNodeList);

  void integrateAndAssembleQueuedMatrixBlocksPart2(
      std::vector<node_t> &admissibleLeafNodeList,
      std::vector<node_t> &queuedLeafNodeList,
      ParallelDataContainer &hMatrixData,
      size_t &finishedAdmissibleLeafNodeCount);

  void checkAcaConvergence(const node_t &blockClusterTreeNode,
                           const size_t index,
                           ParallelDataContainer &hMatrixData,
                           size_t &finishedAdmissibleLeafNodeCount);

  void setAndAnalyseAcaStatus(const size_t index,
                              size_t &finishedAdmissibleLeafNodeCount);

  void analyseAcaStatus(const size_t index,
                        size_t &finishedAdmissibleLeafNodeCount);

  bool addMatrixBlockIntegralsPart1(const node_t &blockClusterTreeNode,
                                    const size_t index);

  bool addMatrixBlockIntegralsPart2(const node_t &blockClusterTreeNode,
                                    const size_t index,
                                    ParallelDataContainer &hMatrixData);

  bool selectMinPivot(const Matrix<ResultType> &vec,
                      const std::vector<std::size_t> &approximationCount,
                      std::size_t &pivot, ModeType modus) const;

  double updateLowRankBlocksAndNorm(const Matrix<ResultType> &newRow,
                                    const Matrix<ResultType> &newCol,
                                    Matrix<ResultType> &A, Matrix<ResultType> &B,
                                    const double currentBlockNorm) const;

  bool checkConvergence(const Matrix<ResultType> &row,
                        const Matrix<ResultType> &col, const double tol) const;

  const Bempp::Space<BasisFunctionType> &m_testSpace;
  const Bempp::Space<BasisFunctionType> &m_trialSpace;
  const shared_ptr<DefaultBlockClusterTreeType> m_blockClusterTree;
  const std::vector<LocalAssembler*> &m_assemblers;
  const std::vector<const DiscreteLinearOp*> &m_sparseTermsToAdd;
  const std::vector<ResultType> &m_denseTermsMultipliers;
  const std::vector<ResultType> &m_sparseTermsMultipliers;
  const double m_eps;
  const double m_zeroTol;
  const unsigned int m_maxRank;
  HMatrixDenseCompressor<ResultType, N> m_hMatrixDenseCompressor;
  shared_ptr<Bempp::LocalDofListsCache<BasisFunctionType>> m_testDofListsCache,
      m_trialDofListsCache;
  mutable tbb::atomic<size_t> m_accessedEntryCount;
  const unsigned int m_generalLocalTestDofCount;

  std::vector<shared_ptr<Fiber::CudaIntegrator<
      BasisFunctionType, KernelType, ResultType,
      CudaBasisFunctionType, CudaKernelType, CudaResultType>>> m_cudaIntegrators;
  std::vector<CudaResultType*> m_results;

  size_t m_maxActiveElemPairCount, m_queuedElemPairCount;
  std::vector<int> m_queuedTestElemPairIndices, m_queuedTrialElemPairIndices;
  std::vector<size_t> m_queuedNodalElemPairCounts;
  size_t m_queuedNodeCount;
  unsigned int m_nextDevice;

  // Further vectors to keep track of the compression state of the
  // admissible leaf nodes
  std::vector<std::vector<size_t>> m_rowApproxCounters;
  std::vector<std::vector<size_t>> m_colApproxCounters;
  std::vector<size_t> m_maxIterations;
  std::vector<Matrix<ResultType>> m_origRows;
  std::vector<Matrix<ResultType>> m_origCols;
  std::vector<double> m_blockNorms;
  std::vector<bool> m_finished;
  std::vector<AcaAlgorithmStateType> m_acaAlgorithmStates;
  std::vector<AcaStatusType> m_acaStatuses;
  std::vector<size_t> m_nextPivots;
  std::vector<ModeType> m_modes;
  std::vector<bool> m_rowTrialPassed;
  std::vector<bool> m_columnTrialPassed;
  std::vector<size_t> m_rowTrialCounts;
  std::vector<size_t> m_columnTrialCounts;
  std::vector<size_t> m_MAX_ROW_TRIAL_COUNTS;
  std::vector<size_t> m_MAX_COLUMN_TRIAL_COUNTS;
  std::vector<bool> m_rowModeZeroTermination;
  std::vector<bool> m_columnModeZeroTermination;

  std::vector<bool> m_acaConverged;
  std::vector<size_t> m_acaNextPivots;
  std::vector<size_t> m_acaIterationCounts;
  std::vector<CrossStatusType> m_acaCrossStatuses;
  std::vector<Matrix<ResultType>> m_rows;
  std::vector<Matrix<ResultType>> m_cols;
  /** \endcond */
};
}

#include "cuda_hmatrix_aca_compressor_impl.hpp"

#endif
