// vi: set et ts=4 sw=2 sts=2:

#ifndef HMAT_HMATRIX_IMPL_HPP
#define HMAT_HMATRIX_IMPL_HPP

#include "hmatrix.hpp"
#include "hmatrix_data.hpp"
#include "hmatrix_dense_data.hpp"
#include "hmatrix_low_rank_data.hpp"
#include "math_helper.hpp"
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/spin_mutex.h>
#include <tbb/task_group.h>

#include <algorithm>
#include <chrono>
#include <cuda_profiler_api.h>

namespace hmat {

template <typename ValueType, int N>
HMatrix<ValueType, N>::HMatrix(
    const shared_ptr<BlockClusterTree<N>> &blockClusterTree, MPI_Comm comm)
    : m_blockClusterTree(blockClusterTree), m_numberOfDenseBlocks(0),
      m_numberOfLowRankBlocks(0), m_memSizeKb(0.0), m_comm(comm) {

  MPI_Comm_size(comm, &m_nproc);
  MPI_Comm_rank(comm, &m_rank);
}

template <typename ValueType, int N>
HMatrix<ValueType, N>::HMatrix(
    const shared_ptr<BlockClusterTree<N>> &blockClusterTree,
    const HMatrixCompressor<ValueType, N> &hMatrixCompressor, MPI_Comm comm)
    : HMatrix<ValueType, N>(blockClusterTree, comm) {
  initialize(hMatrixCompressor);
}

template <typename ValueType, int N>
HMatrix<ValueType, N>::HMatrix(
  ParallelDataContainer &hMatrixData,
  const shared_ptr<BlockClusterTree<N>> &blockClusterTree,
  bool coarsening, double coarsening_accuracy)
  : m_hMatrixData(hMatrixData), m_blockClusterTree(blockClusterTree),
    m_numberOfDenseBlocks(0), m_numberOfLowRankBlocks(0), m_memSizeKb(0.0) {

  typedef decltype(m_blockClusterTree->root()) node_t;

  std::function<void(const node_t &node)> coarsenFun =
      [&](const node_t &node) {
        if (!node->isLeaf()) {
          tbb::task_group g;
          g.run([&] { coarsenFun(node->child(0)); });
          g.run([&] { coarsenFun(node->child(1)); });
          g.run([&] { coarsenFun(node->child(2)); });
          g.run_and_wait([&] { coarsenFun(node->child(3)); });

          // Now do a coarsen step
          coarsen_impl(node, coarsening_accuracy);
        }
      };

  // Start the coarsening
  if (coarsening)
    coarsenFun(m_blockClusterTree->root());

  // Compute statistics
  for (auto &elem : m_hMatrixData) {
    if (!elem.second)
      continue;
    if (elem.second->type() == DENSE)
      m_numberOfDenseBlocks++;
    else
      m_numberOfLowRankBlocks++;
    m_memSizeKb += elem.second->memSizeKb();
  }
}

template <typename ValueType, int N>
std::size_t HMatrix<ValueType, N>::rows() const {
  return m_blockClusterTree->rows();
}

template <typename ValueType, int N>
std::size_t HMatrix<ValueType, N>::columns() const {
  return m_blockClusterTree->columns();
}

template <typename ValueType, int N>
double HMatrix<ValueType, N>::frobeniusNorm() const {
  return frobeniusNorm_impl(m_blockClusterTree->root());
}

template <typename ValueType, int N>
int HMatrix<ValueType, N>::numberOfDenseBlocks() const {
  return m_numberOfDenseBlocks;
}

template <typename ValueType, int N>
int HMatrix<ValueType, N>::numberOfLowRankBlocks() const {
  return m_numberOfLowRankBlocks;
}

template <typename ValueType, int N>
int HMatrix<ValueType, N>::numberOfBlocks() const {
  return m_numberOfLowRankBlocks + m_numberOfDenseBlocks;
}

template <typename ValueType, int N>
double HMatrix<ValueType, N>::memSizeKb() const {
  return m_memSizeKb;
}

template <typename ValueType, int N>
void HMatrix<ValueType, N>::initialize(
    const HMatrixCompressor<ValueType, N> &hMatrixCompressor) {

  reset();

  auto leafNodes = this->m_blockClusterTree->leafNodes();

  // Type of a block cluster tree node.
  typedef decltype(m_blockClusterTree->root()) node_t;

  // Sorting by size not necessary any more as a natural size
  // ordering is induced by the breadth-first search of the
  // leafs function that returns the leaf vector.
  /*

    // Sort the leaf nodes by size.
    std::stable_sort(
        leafNodes.begin(), leafNodes.end(),
        [&](const node_t &a, const node_t &b) -> bool {
          auto rowRange_a = a->data().rowClusterTreeNode->data().indexRange;
          auto colRange_a = a->data().columnClusterTreeNode->data().indexRange;
          auto rowRange_b = b->data().rowClusterTreeNode->data().indexRange;
          auto colRange_b = b->data().columnClusterTreeNode->data().indexRange;

          std::size_t size_a = std::max(rowRange_a[1] - rowRange_a[0],
                                        colRange_a[1] - colRange_a[0]);
          std::size_t size_b = std::max(rowRange_b[1] - rowRange_b[0],
                                        colRange_b[1] - colRange_b[0]);

          // Sort in descending order
          return size_a > size_b;

        });
  */

  for (int index = m_rank; index < leafNodes.size(); index += m_nproc)
    m_myLeafs.push_back(leafNodes[index]);

  std::size_t numberOfLeafs = m_myLeafs.size();

  cudaProfilerStart();

  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, numberOfLeafs),
                    [&](const tbb::blocked_range<std::size_t> &r) {
                      for (auto index = r.begin(); index != r.end(); index++) {
                        shared_ptr<HMatrixData<ValueType>> nodeData;
                        hMatrixCompressor.compressBlock(*m_myLeafs[index],
                                                        nodeData);
                        m_hMatrixData[m_myLeafs[index]] = nodeData;
                      }
                    });

  cudaProfilerStop();

  // Compute statistics
  for (auto &elem : m_hMatrixData) {
    if (!elem.second)
      continue;
    if (elem.second->type() == DENSE)
      m_numberOfDenseBlocks++;
    else
      m_numberOfLowRankBlocks++;
    m_memSizeKb += elem.second->memSizeKb();
  }
  MPI_Barrier(m_comm);
}

template <typename ValueType, int N> void HMatrix<ValueType, N>::reset() {
  m_hMatrixData.clear();
}

template <typename ValueType, int N>
bool HMatrix<ValueType, N>::isInitialized() const {
  return (!m_hMatrixData.empty());
}

template <typename ValueType, int N>
shared_ptr<const BlockClusterTree<N>>
HMatrix<ValueType, N>::blockClusterTree() const {
  return this->m_blockClusterTree;
}

template <typename ValueType, int N>
shared_ptr<const HMatrixData<ValueType>> HMatrix<ValueType, N>::data(
    shared_ptr<const BlockClusterTreeNode<N>> node) const {
  return this->m_hMatrixData.at(
      const_pointer_cast<BlockClusterTreeNode<N>>(node));
}

template <typename ValueType, int N>
Matrix<ValueType> HMatrix<ValueType, N>::permuteMatToHMatDofs(
    const Eigen::Ref<Matrix<ValueType>> &mat,
    RowColSelector rowOrColumn) const {

  Matrix<ValueType> permutedDofs(mat.rows(), mat.cols());

  shared_ptr<const ClusterTree<N>> clusterTree;

  if (rowOrColumn == ROW)
    clusterTree = m_blockClusterTree->rowClusterTree();
  else
    clusterTree = m_blockClusterTree->columnClusterTree();

  for (std::size_t i = 0; i < mat.rows(); ++i) {
    auto permutedIndex = clusterTree->mapOriginalDofToHMatDof(i);
    for (std::size_t j = 0; j < mat.cols(); ++j)
      permutedDofs(permutedIndex, j) = mat(i, j);
  }

  return permutedDofs;
}

template <typename ValueType, int N>
Matrix<ValueType> HMatrix<ValueType, N>::permuteMatToOriginalDofs(
    const Eigen::Ref<Matrix<ValueType>> &mat,
    RowColSelector rowOrColumn) const {

  Matrix<ValueType> originalDofs(mat.rows(), mat.cols());

  shared_ptr<const ClusterTree<N>> clusterTree;

  if (rowOrColumn == ROW)
    clusterTree = m_blockClusterTree->rowClusterTree();
  else
    clusterTree = m_blockClusterTree->columnClusterTree();

  for (std::size_t i = 0; i < mat.rows(); ++i) {
    auto originalIndex = clusterTree->mapHMatDofToOriginalDof(i);
    for (std::size_t j = 0; j < mat.cols(); ++j)
      originalDofs(originalIndex, j) = mat(i, j);
  }

  return originalDofs;
}

template <typename ValueType, int N>
void HMatrix<ValueType, N>::apply(const Eigen::Ref<Matrix<ValueType>> &X,
                                  Eigen::Ref<Matrix<ValueType>> Y,
                                  TransposeMode trans, ValueType alpha,
                                  ValueType beta) const {

  if (beta == ValueType(0) || m_rank != 0)
    Y.setZero();
  else {
    Y *= beta;
  }

  Matrix<ValueType> xPermuted;
  Matrix<ValueType> yPermuted;

  if (trans == TransposeMode::NOTRANS || trans == TransposeMode::CONJ) {

    xPermuted = alpha * permuteMatToHMatDofs(X, COL);
    yPermuted = permuteMatToHMatDofs(Y, ROW);

  } else {

    xPermuted = alpha * permuteMatToHMatDofs(X, ROW);
    yPermuted = permuteMatToHMatDofs(Y, COL);
  }

  std::size_t numberOfLeafs = m_myLeafs.size();
  auto cols = xPermuted.cols();

  // Hack because Eigen ref does not like const.
  // Eigen::Ref<Matrix<ValueType>> x_no_const =
  // const_cast<Eigen::Ref<Matrix<ValueType>> &>(xPermuted);

  // Initialize one mutex for each output row
  Matrix<tbb::spin_mutex> mutex(Y.rows(), 1);

  tbb::parallel_for(
      tbb::blocked_range<std::size_t>(0, numberOfLeafs),
      [&](const tbb::blocked_range<std::size_t> &r) {
        for (auto index = r.begin(); index != r.end(); ++index) {
          auto rowRange =
              m_myLeafs[index]->data().rowClusterTreeNode->data().indexRange;
          auto colRange =
              m_myLeafs[index]->data().columnClusterTreeNode->data().indexRange;

          if (trans == TransposeMode::NOTRANS || trans == TransposeMode::CONJ) {
            Matrix<ValueType> x_data = xPermuted.block(
                colRange[0], 0, colRange[1] - colRange[0], cols);
            Matrix<ValueType> result =
                Matrix<ValueType>::Zero(rowRange[1] - rowRange[0], cols);
            this->m_hMatrixData.at(m_myLeafs[index])
                ->apply(x_data, result, trans, 1, 1);

            for (int i = rowRange[0]; i < rowRange[1]; ++i) {
              tbb::spin_mutex::scoped_lock lock(mutex(i, 0));
              for (int j = 0; j < cols; ++j)
                yPermuted(i, j) += result(i - rowRange[0], j);
            }
          } else {
            Eigen::Ref<Matrix<ValueType>> x_data = xPermuted.block(
                rowRange[0], 0, rowRange[1] - rowRange[0], cols);
            Matrix<ValueType> result =
                Matrix<ValueType>::Zero(colRange[1] - colRange[0], cols);
            this->m_hMatrixData.at(m_myLeafs[index])
                ->apply(x_data, result, trans, 1, 1);

            for (int i = colRange[0]; i < colRange[1]; ++i) {
              tbb::spin_mutex::scoped_lock lock(mutex(i, 0));
              for (int j = 0; j < cols; ++j)
                yPermuted(i, j) += result(i - colRange[0], j);
            }
          }
        }
      });

  if (trans == TransposeMode::NOTRANS || trans == TransposeMode::CONJ)
    Y = this->permuteMatToOriginalDofs(yPermuted, ROW);
  else
    Y = this->permuteMatToOriginalDofs(yPermuted, COL);

  if (m_nproc > 1) {
    // Perform an all reduce operation if there is more than one proc
    Matrix<ValueType> globalY(Y.rows(), Y.cols());
    MPI_Allreduce(Y.data(), globalY.data(), Y.rows() * Y.cols(),
                  MpiTrait<ValueType>().type, MPI_SUM, m_comm);
    Y = globalY;
  }
}

template <typename ValueType, int N>
double HMatrix<ValueType, N>::frobeniusNorm_impl(
    const shared_ptr<BlockClusterTreeNode<N>> &node) const {

  if (node->isLeaf())
    return m_hMatrixData.at(node)->frobeniusNorm();

  tbb::task_group g;

  double result = 0;

  double res0;
  double res1;
  double res2;
  double res3;

  g.run([&] { res0 = frobeniusNorm_impl(node->child(0)); });
  g.run([&] { res1 = frobeniusNorm_impl(node->child(1)); });
  g.run([&] { res2 = frobeniusNorm_impl(node->child(2)); });
  g.run_and_wait([&] { res3 = frobeniusNorm_impl(node->child(3)); });

  result = std::sqrt(res0 * res0 + res1 * res1 + res2 * res2 + res3 * res3);

  return result;
}
}

#endif
