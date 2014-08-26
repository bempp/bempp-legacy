// Copyright (C) 2011-2012 by the BEM++ Authors
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

#include "bempp/common/config_ahmed.hpp"
#include "bempp/common/config_trilinos.hpp"
#ifdef WITH_AHMED

#include "discrete_aca_boundary_operator.hpp"
#include "../common/shared_ptr.hpp"

#include "ahmed_aux.hpp"
#include "aca_approximate_lu_inverse.hpp"

#include "../common/chunk_statistics.hpp"
#include "../common/complex_aux.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/serial_blas_region.hpp"

#include <fstream>
#include <iostream>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/type_traits/is_complex.hpp>

#include <tbb/blocked_range.h>
#include <tbb/concurrent_queue.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>

#ifdef WITH_TRILINOS
#include <Thyra_DefaultSpmdVectorSpace_decl.hpp>
#endif

namespace Bempp {

void dumpblcluster(const blcluster *bl, const std::string &indent) {
  std::cout << indent << bl << " " << bl->getb1() << " " << bl->getb2() << " "
            << bl->getn1() << " " << bl->getn2() << "\n";
  std::string sonIndent = indent + "  ";
  for (int r = 0; r < bl->getnrs(); ++r)
    for (int c = 0; c < bl->getncs(); ++c) {
      blcluster *son = bl->getson(r, c);
      if (son)
        dumpblcluster(son, sonIndent);
      else
        std::cout << sonIndent << "0\n";
    }
}

namespace {

template <typename ValueType> class MblockMultiplicationLoopBody {
  typedef mblock<typename AhmedTypeTraits<ValueType>::Type> AhmedMblock;

public:
  typedef tbb::concurrent_queue<size_t> LeafClusterIndexQueue;

  MblockMultiplicationLoopBody(TranspositionMode trans, ValueType multiplier,
                               arma::Col<ValueType> &x, arma::Col<ValueType> &y,
                               AhmedLeafClusterArray &leafClusters,
                               boost::shared_array<AhmedMblock *> blocks,
                               LeafClusterIndexQueue &leafClusterIndexQueue,
                               std::vector<ChunkStatistics> &stats)
      : m_trans(trans), m_multiplier(multiplier), m_x(x), m_local_y(y),
        m_leafClusters(leafClusters), m_blocks(blocks),
        m_leafClusterIndexQueue(leafClusterIndexQueue), m_stats(stats) {
    if (trans != NO_TRANSPOSE && trans != TRANSPOSE &&
        trans != CONJUGATE_TRANSPOSE)
      throw std::invalid_argument(
          "MblockMultiplicationLoopBody::MblockMultiplicationLoopBody(): "
          "unsupported transposition mode");
    // m_local_y.fill(static_cast<ValueType>(0.));
  }

  MblockMultiplicationLoopBody(MblockMultiplicationLoopBody &other, tbb::split)
      : m_trans(other.m_trans), m_multiplier(other.m_multiplier),
        m_x(other.m_x), m_local_y(other.m_local_y.n_rows),
        m_leafClusters(other.m_leafClusters), m_blocks(other.m_blocks),
        m_leafClusterIndexQueue(other.m_leafClusterIndexQueue),
        m_stats(other.m_stats) {
    m_local_y.fill(static_cast<ValueType>(0.));
  }

  template <typename Range> void operator()(const Range &r) {
    for (typename Range::const_iterator i = r.begin(); i != r.end(); ++i) {
      size_t leafClusterIndex = -1;
      if (!m_leafClusterIndexQueue.try_pop(leafClusterIndex)) {
        std::cerr << "MblockMultiplicationLoopBody::operator(): "
                     "Warning: try_pop failed; this shouldn't happen!"
                  << std::endl;
        continue;
      }
      m_stats[leafClusterIndex].valid = true;
      m_stats[leafClusterIndex].chunkStart = r.begin();
      m_stats[leafClusterIndex].chunkSize = r.size();
      m_stats[leafClusterIndex].startTime = tbb::tick_count::now();

      blcluster *cluster = m_leafClusters[leafClusterIndex];
      if (m_trans == NO_TRANSPOSE)
        m_blocks[cluster->getidx()]->mltaVec(
            ahmedCast(m_multiplier), ahmedCast(&m_x(cluster->getb2())),
            ahmedCast(&m_local_y(cluster->getb1())));
      else if (m_trans == TRANSPOSE)
        m_blocks[cluster->getidx()]->mltatVec(
            ahmedCast(m_multiplier), ahmedCast(&m_x(cluster->getb1())),
            ahmedCast(&m_local_y(cluster->getb2())));
      else // m_trans == CONJUGATE_TRANSPOSE)
        m_blocks[cluster->getidx()]->mltahVec(
            ahmedCast(m_multiplier), ahmedCast(&m_x(cluster->getb1())),
            ahmedCast(&m_local_y(cluster->getb2())));
      m_stats[leafClusterIndex].endTime = tbb::tick_count::now();
    }
  }

  void join(const MblockMultiplicationLoopBody &other) {
    m_local_y += other.m_local_y;
  }

private:
  TranspositionMode m_trans;
  ValueType m_multiplier;
  arma::Col<ValueType> &m_x;

public:
  arma::Col<ValueType> m_local_y;

private:
  AhmedLeafClusterArray &m_leafClusters;
  boost::shared_array<AhmedMblock *> m_blocks;
  LeafClusterIndexQueue &m_leafClusterIndexQueue;
  std::vector<ChunkStatistics> &m_stats;
};

bool areEqual(const blcluster *op1, const blcluster *op2) {
  if (!op1 || !op2)
    return (!op1 && !op2);
  if (op1->getnrs() != op2->getnrs() || op1->getncs() != op2->getncs())
    return false;
  if (op1->isleaf())
    return (op1->getb1() == op2->getb1() && op1->getb2() == op2->getb2() &&
            op1->getn1() == op2->getn1() && op1->getn2() == op2->getn2());
  for (unsigned int rs = 0; rs < op1->getnrs(); ++rs)
    for (unsigned int cs = 0; cs < op1->getncs(); ++cs)
      if (!areEqual(op1->getson(rs, cs), op2->getson(rs, cs)))
        return false;
  return true;
}

inline void mltaSyHVec(double d, blcluster *bl, mblock<double> **A, double *x,
                       double *y) {
  mltaHeHVec(d, bl, A, x, y);
}

inline void mltaSyHVec(float d, blcluster *bl, mblock<float> **A, float *x,
                       float *y) {
  mltaHeHVec(d, bl, A, x, y);
}

inline void mltaSyHVec(scomp d, blcluster *bl, mblock<scomp> **A, scomp *x,
                       scomp *y) {
  throw std::runtime_error("mltaSyHVec(): the overload of this function for "
                           "complex single-precision numbers does not exist "
                           "in AHMED");
}

// For real types, SyHh = SyH
inline void mltaSyHhVec(double d, blcluster *bl, mblock<double> **A, double *x,
                        double *y) {
  mltaSyHVec(d, bl, A, x, y);
}

inline void mltaSyHhVec(float d, blcluster *bl, mblock<float> **A, float *x,
                        float *y) {
  mltaSyHVec(d, bl, A, x, y);
}

inline void mltaSyHhVec(scomp d, blcluster *bl, mblock<scomp> **A, scomp *x,
                        scomp *y) {
  throw std::runtime_error("mltaSyHhVec(): the overload of this function for "
                           "complex single-precision numbers does not exist "
                           "in AHMED");
}

} // namespace

template <typename ValueType>
DiscreteAcaBoundaryOperator<ValueType>::DiscreteAcaBoundaryOperator(
    unsigned int rowCount, unsigned int columnCount, double eps_,
    int maximumRank_, int symmetry_,
    const shared_ptr<const AhmedBemBlcluster> &blockCluster_,
    const AhmedMblockArray &blocks_, const IndexPermutation &domainPermutation_,
    const IndexPermutation &rangePermutation_,
    const ParallelizationOptions &parallelizationOptions_,
    const std::vector<AhmedConstMblockArray> &sharedBlocks_)
    :
#ifdef WITH_TRILINOS
      m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(columnCount)),
      m_rangeSpace(Thyra::defaultSpmdVectorSpace<ValueType>(rowCount)),
#else
      m_rowCount(rowCount),
      m_columnCount(columnCount),
#endif
      m_eps(eps_), m_maximumRank(maximumRank_), m_symmetry(symmetry_),
      m_blockCluster(blockCluster_), m_blocks(blocks_),
      m_domainPermutation(domainPermutation_),
      m_rangePermutation(rangePermutation_),
      m_parallelizationOptions(parallelizationOptions_),
      m_sharedBlocks(sharedBlocks_) {
  if (eps_ <= 0 || eps_ > 1)
    std::cout << "DiscreteAcaBoundaryOperator::DiscreteAcaBoundaryOperator(): "
                 "warning: suspicious value of eps (" << eps_ << ")"
              << std::endl;
  if (maximumRank_ < 0)
    std::cout << "DiscreteAcaBoundaryOperator::DiscreteAcaBoundaryOperator(): "
                 "warning: suspicious value of maximumRank (" << maximumRank_
              << ")" << std::endl;
}

template <typename ValueType>
DiscreteAcaBoundaryOperator<ValueType>::DiscreteAcaBoundaryOperator(
    unsigned int rowCount, unsigned int columnCount, int maximumRank_,
    int symmetry_, std::unique_ptr<const AhmedBemBlcluster> blockCluster_,
    AhmedMblockArray blocks_, const IndexPermutation &domainPermutation_,
    const IndexPermutation &rangePermutation_,
    const ParallelizationOptions &parallelizationOptions_,
    const std::vector<AhmedConstMblockArray> &sharedBlocks_)
    :
#ifdef WITH_TRILINOS
      m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(columnCount)),
      m_rangeSpace(Thyra::defaultSpmdVectorSpace<ValueType>(rowCount)),
#else
      m_rowCount(rowCount),
      m_columnCount(columnCount),
#endif
      m_eps(1e-4), // just a guess...
      m_maximumRank(maximumRank_), m_symmetry(symmetry_),
      m_blockCluster(blockCluster_.release()), m_blocks(blocks_),
      m_domainPermutation(domainPermutation_),
      m_rangePermutation(rangePermutation_),
      m_parallelizationOptions(parallelizationOptions_),
      m_sharedBlocks(sharedBlocks_) {
}

template <typename ValueType>
DiscreteAcaBoundaryOperator<ValueType>::DiscreteAcaBoundaryOperator(
    unsigned int rowCount, unsigned int columnCount, int maximumRank_,
    int symmetry_, const shared_ptr<const AhmedBemBlcluster> &blockCluster_,
    const AhmedMblockArray &blocks_, const IndexPermutation &domainPermutation_,
    const IndexPermutation &rangePermutation_,
    const ParallelizationOptions &parallelizationOptions_,
    const std::vector<AhmedConstMblockArray> &sharedBlocks_)
    :
#ifdef WITH_TRILINOS
      m_domainSpace(Thyra::defaultSpmdVectorSpace<ValueType>(columnCount)),
      m_rangeSpace(Thyra::defaultSpmdVectorSpace<ValueType>(rowCount)),
#else
      m_rowCount(rowCount),
      m_columnCount(columnCount),
#endif
      m_eps(1e-4), // just a guess...
      m_maximumRank(maximumRank_), m_symmetry(symmetry_),
      m_blockCluster(blockCluster_), m_blocks(blocks_),
      m_domainPermutation(domainPermutation_),
      m_rangePermutation(rangePermutation_),
      m_parallelizationOptions(parallelizationOptions_),
      m_sharedBlocks(sharedBlocks_) {
}

template <typename ValueType>
arma::Mat<ValueType> DiscreteAcaBoundaryOperator<ValueType>::asMatrix() const {
  const unsigned int nRows = rowCount();
  const unsigned int nCols = columnCount();
  const blcluster *blockCluster = m_blockCluster.get();
  blcluster *nonconstBlockCluster = const_cast<blcluster *>(blockCluster);

  arma::Mat<ValueType> permutedOutput(nRows, nCols);
  permutedOutput.fill(0.);
  arma::Col<ValueType> unit(nCols);
  unit.fill(0.);

  for (unsigned int col = 0; col < nCols; ++col) {
    if (col > 0)
      unit(col - 1) = 0.;
    unit(col) = 1.;
    if (m_symmetry & SYMMETRIC)
      mltaSyHVec(1., nonconstBlockCluster, m_blocks.get(),
                 ahmedCast(unit.memptr()),
                 ahmedCast(permutedOutput.colptr(col)));
    else if (m_symmetry & HERMITIAN)
      mltaHeHVec(1., nonconstBlockCluster, m_blocks.get(),
                 ahmedCast(unit.memptr()),
                 ahmedCast(permutedOutput.colptr(col)));
    else
      mltaGeHVec(1., nonconstBlockCluster, m_blocks.get(),
                 ahmedCast(unit.memptr()),
                 ahmedCast(permutedOutput.colptr(col)));
  }

  arma::Mat<ValueType> output(nRows, nCols);
  for (unsigned int col = 0; col < nCols; ++col)
    for (unsigned int row = 0; row < nRows; ++row)
      output(row, col) = permutedOutput(m_rangePermutation.permuted(row),
                                        m_domainPermutation.permuted(col));

  return output;
}

template <typename ValueType>
unsigned int DiscreteAcaBoundaryOperator<ValueType>::rowCount() const {
#ifdef WITH_TRILINOS
  return m_rangeSpace->dim();
#else
  return m_rowCount;
#endif
}

template <typename ValueType>
unsigned int DiscreteAcaBoundaryOperator<ValueType>::columnCount() const {
#ifdef WITH_TRILINOS
  return m_domainSpace->dim();
#else
  return m_columnCount;
#endif
}

template <typename ValueType>
void DiscreteAcaBoundaryOperator<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, arma::Mat<ValueType> &block) const {
  throw std::runtime_error("DiscreteAcaBoundaryOperator::"
                           "addBlock(): not implemented yet");
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>>
DiscreteAcaBoundaryOperator<ValueType>::asDiscreteAcaBoundaryOperator(
    double eps, int maximumRank, bool interleave) const {
  return this->shared_from_this(); // this-> needed for template name
                                   // resolution.
}

template <typename ValueType>
const DiscreteAcaBoundaryOperator<ValueType> &
DiscreteAcaBoundaryOperator<ValueType>::castToAca(
    const DiscreteBoundaryOperator<ValueType> &discreteOperator) {
  return dynamic_cast<const DiscreteAcaBoundaryOperator<ValueType> &>(
      discreteOperator);
}

template <typename ValueType>
shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>>
DiscreteAcaBoundaryOperator<ValueType>::castToAca(const shared_ptr<
    const DiscreteBoundaryOperator<ValueType>> &discreteOperator) {
  shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>> result =
      boost::dynamic_pointer_cast<const DiscreteAcaBoundaryOperator<ValueType>>(
          discreteOperator);
  if (result.get() == 0 && (discreteOperator.get() != 0))
    throw std::bad_cast();
  return result;
}

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
DiscreteAcaBoundaryOperator<ValueType>::domain() const {
  return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
DiscreteAcaBoundaryOperator<ValueType>::range() const {
  return m_rangeSpace;
}

template <typename ValueType>
bool DiscreteAcaBoundaryOperator<ValueType>::opSupportedImpl(
    Thyra::EOpTransp M_trans) const {
  // TODO: implement remaining variant (conjugate)
  return (M_trans == Thyra::NOTRANS || M_trans == Thyra::TRANS ||
          M_trans == Thyra::CONJTRANS);
}
#endif // WITH_TRILINOS

template <typename ValueType>
void DiscreteAcaBoundaryOperator<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const arma::Col<ValueType> &x_in,
    arma::Col<ValueType> &y_inout, const ValueType alpha,
    const ValueType beta) const {
  if (trans != NO_TRANSPOSE && trans != TRANSPOSE &&
      trans != CONJUGATE_TRANSPOSE)
    throw std::runtime_error(
        "DiscreteAcaBoundaryOperator::applyBuiltInImpl(): "
        "transposition modes other than NO_TRANSPOSE, TRANSPOSE and "
        "CONJUGATE_TRANSPOSE are not supported");
  bool transposed = (trans & TRANSPOSE);

  const blcluster *blockCluster = m_blockCluster.get();
  blcluster *nonconstBlockCluster = const_cast<blcluster *>(blockCluster);

  if ((!transposed &&
       (columnCount() != x_in.n_rows || rowCount() != y_inout.n_rows)) ||
      (transposed &&
       (rowCount() != x_in.n_rows || columnCount() != y_inout.n_rows)))
    throw std::invalid_argument(
        "DiscreteAcaBoundaryOperator::applyBuiltInImpl(): "
        "incorrect vector length");

  if (beta == static_cast<ValueType>(0.))
    y_inout.fill(static_cast<ValueType>(0.));
  else
    y_inout *= beta;

  arma::Col<ValueType> permutedArgument;
  if (!transposed)
    m_domainPermutation.permuteVector(x_in, permutedArgument);
  else
    m_rangePermutation.permuteVector(x_in, permutedArgument);

  arma::Col<ValueType> permutedResult;
  if (!transposed)
    m_rangePermutation.permuteVector(y_inout, permutedResult);
  else
    m_domainPermutation.permuteVector(y_inout, permutedResult);

  if (m_symmetry & SYMMETRIC) {
    // TODO: parallelise
    if (trans == NO_TRANSPOSE || trans == TRANSPOSE)
      mltaSyHVec(ahmedCast(alpha), nonconstBlockCluster, m_blocks.get(),
                 ahmedCast(permutedArgument.memptr()),
                 ahmedCast(permutedResult.memptr()));
    else // trans == CONJUGATE_TRANSPOSE
      mltaSyHhVec(ahmedCast(alpha), nonconstBlockCluster, m_blocks.get(),
                  ahmedCast(permutedArgument.memptr()),
                  ahmedCast(permutedResult.memptr()));
  } else if (m_symmetry & HERMITIAN) {
    // NO_TRANSPOSE and CONJUGATE_TRANSPOSE are equivalent
    if (trans == NO_TRANSPOSE || trans == CONJUGATE_TRANSPOSE)
      mltaHeHVec(ahmedCast(alpha), nonconstBlockCluster, m_blocks.get(),
                 ahmedCast(permutedArgument.memptr()),
                 ahmedCast(permutedResult.memptr()));
    else { // trans == TRANSPOSE
      // alpha A^T x + beta y = (alpha^* A^H x^* + beta^* y^*)^*
      // = (alpha^* A x^* + beta^* y^*)^*
      permutedArgument = arma::conj(permutedArgument);
      permutedResult = arma::conj(permutedResult);
      ValueType alphaConj = conj(alpha);
      mltaHeHVec(ahmedCast(alphaConj), nonconstBlockCluster, m_blocks.get(),
                 ahmedCast(permutedArgument.memptr()),
                 ahmedCast(permutedResult.memptr()));
      permutedResult = arma::conj(permutedResult);
    }
  } else {
    //        if (trans == NO_TRANSPOSE)
    //            mltaGeHVec(ahmedCast(alpha), nonconstBlockCluster,
    // m_blocks.get(),
    //                       ahmedCast(permutedArgument.memptr()),
    //                       ahmedCast(permutedResult.memptr()));
    //        else // trans == CONJUGATE_TRANSPOSE
    //            mltaGeHhVec(ahmedCast(alpha), nonconstBlockCluster,
    // m_blocks.get(),
    //                        ahmedCast(permutedArgument.memptr()),
    //                        ahmedCast(permutedResult.memptr()));

    AhmedLeafClusterArray leafClusters(nonconstBlockCluster);
    leafClusters.sortAccordingToClusterSize();
    const size_t leafClusterCount = leafClusters.size();

    int maxThreadCount = 1;
    if (!m_parallelizationOptions.isOpenClEnabled()) {
      if (m_parallelizationOptions.maxThreadCount() ==
          ParallelizationOptions::AUTO)
        maxThreadCount = tbb::task_scheduler_init::automatic;
      else
        maxThreadCount = m_parallelizationOptions.maxThreadCount();
    }
    tbb::task_scheduler_init scheduler(maxThreadCount);

    std::vector<ChunkStatistics> chunkStats(leafClusterCount);

    typedef MblockMultiplicationLoopBody<ValueType> Body;
    typename Body::LeafClusterIndexQueue leafClusterIndexQueue;
    for (size_t i = 0; i < leafClusterCount; ++i)
      leafClusterIndexQueue.push(i);

    // std::cout << "----------------------------\nperm Arg\n" <<
    // permutedArgument;
    // std::cout << "perm Res\n" << permutedResult;
    Body body(trans, alpha, permutedArgument, permutedResult, leafClusters,
              m_blocks, leafClusterIndexQueue, chunkStats);
    {
      Fiber::SerialBlasRegion region;
      tbb::parallel_reduce(tbb::blocked_range<size_t>(0, leafClusterCount),
                           body);
    }
    permutedResult = body.m_local_y;
  }
  if (!transposed)
    m_rangePermutation.unpermuteVector(permutedResult, y_inout);
  else
    m_domainPermutation.unpermuteVector(permutedResult, y_inout);
}

template <typename ValueType>
void DiscreteAcaBoundaryOperator<ValueType>::makeAllMblocksDense() {
  for (unsigned int i = 0; i < m_blockCluster->nleaves(); ++i)
    if (m_blocks[i]->isLrM())
      m_blocks[i]->convLrM_toGeM();
}

template <typename ValueType>
double DiscreteAcaBoundaryOperator<ValueType>::eps() const {
  return m_eps;
}

template <typename ValueType>
int DiscreteAcaBoundaryOperator<ValueType>::maximumRank() const {
  return m_maximumRank;
}

template <typename ValueType>
int DiscreteAcaBoundaryOperator<ValueType>::actualMaximumRank() const {
  // const_cast because Ahmed is not const-correct
  return Hmax_rank(const_cast<AhmedBemBlcluster *>(m_blockCluster.get()),
                   m_blocks.get());
}

template <typename ValueType>
int DiscreteAcaBoundaryOperator<ValueType>::symmetry() const {
  return m_symmetry;
}

template <typename ValueType>
shared_ptr<
    const typename DiscreteAcaBoundaryOperator<ValueType>::AhmedBemBlcluster>
DiscreteAcaBoundaryOperator<ValueType>::blockCluster() const {
  return m_blockCluster;
}

template <typename ValueType>
typename DiscreteAcaBoundaryOperator<ValueType>::AhmedMblockArray
DiscreteAcaBoundaryOperator<ValueType>::blocks() const {
  return m_blocks;
}

template <typename ValueType>
size_t DiscreteAcaBoundaryOperator<ValueType>::blockCount() const {
  return m_blockCluster->nleaves();
}

template <typename ValueType>
const IndexPermutation &
DiscreteAcaBoundaryOperator<ValueType>::domainPermutation() const {
  return m_domainPermutation;
}

template <typename ValueType>
const IndexPermutation &
DiscreteAcaBoundaryOperator<ValueType>::rangePermutation() const {
  return m_rangePermutation;
}

template <typename ValueType>
const ParallelizationOptions &
DiscreteAcaBoundaryOperator<ValueType>::parallelizationOptions() const {
  return m_parallelizationOptions;
}

template <typename ValueType>
std::vector<
    typename DiscreteAcaBoundaryOperator<ValueType>::AhmedConstMblockArray>
DiscreteAcaBoundaryOperator<ValueType>::sharedBlocks() const {
  return m_sharedBlocks;
}

// Global routines

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>>
acaOperatorSum(const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op1,
               const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op2,
               double eps, int maximumRank) {
  if (!op1 || !op2)
    throw std::invalid_argument("acaOperatorSum(): "
                                "both operands must be non-null");
  if (eps <= 0)
    throw std::invalid_argument("acaOperatorSum(): eps must be positive");
  if (maximumRank < 1)
    maximumRank = 1;
  shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>> acaOp1 =
      DiscreteAcaBoundaryOperator<ValueType>::castToAca(op1);
  shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>> acaOp2 =
      DiscreteAcaBoundaryOperator<ValueType>::castToAca(op2);
  if (acaOp1->symmetry() != acaOp2->symmetry())
    throw std::runtime_error("acaOperatorSum(): addition of two H-matrices "
                             "of different symmetry is not supported yet");
  if (!areEqual(acaOp1->blockCluster().get(), acaOp2->blockCluster().get()))
    throw std::invalid_argument("acaOperatorSum(): block cluster trees of "
                                "both operands must be identical");
  if (acaOp1->m_domainPermutation != acaOp2->m_domainPermutation ||
      acaOp1->m_rangePermutation != acaOp2->m_rangePermutation)
    throw std::invalid_argument("acaOperatorSum(): domain and range "
                                "index permutations of "
                                "both operands must be identical");

  typedef typename DiscreteAcaBoundaryOperator<ValueType>::AhmedBemBlcluster
  AhmedBemBlcluster;
  typedef typename DiscreteAcaBoundaryOperator<ValueType>::AhmedMblock
  AhmedMblock;

  shared_ptr<const AhmedBemBlcluster> sumBlockCluster = acaOp1->blockCluster();
  blcluster *nonConstSumBlockCluster =
      const_cast<blcluster *>(            // AHMED is not const-correct
          static_cast<const blcluster *>( // safe -- upcast
              sumBlockCluster.get()));
  boost::shared_array<AhmedMblock *> sumBlocks =
      allocateAhmedMblockArray<ValueType>(sumBlockCluster.get());
  copyH(nonConstSumBlockCluster, acaOp1->m_blocks.get(), sumBlocks.get());
  addGeHGeH(nonConstSumBlockCluster, sumBlocks.get(), acaOp2->m_blocks.get(),
            eps, maximumRank);
  shared_ptr<const DiscreteBoundaryOperator<ValueType>> result(
      new DiscreteAcaBoundaryOperator<ValueType>(
          acaOp1->rowCount(), acaOp1->columnCount(), eps, maximumRank,
          acaOp1->m_symmetry & acaOp2->m_symmetry, sumBlockCluster, sumBlocks,
          acaOp1->m_domainPermutation, acaOp2->m_rangePermutation,
          acaOp1->m_parallelizationOptions));
  return result;
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>> scaledAcaOperator(
    const ValueType &multiplier,
    const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op) {
  if (!op)
    throw std::invalid_argument("scaledAcaOperator(): "
                                "operand must not be null");

  typedef typename DiscreteAcaBoundaryOperator<ValueType>::AhmedBemBlcluster
  AhmedBemBlcluster;
  typedef typename DiscreteAcaBoundaryOperator<ValueType>::AhmedMblock
  AhmedMblock;

  typedef typename AhmedTypeTraits<ValueType>::Type AhmedValueType;
  const AhmedValueType ahmedMultiplier = ahmedCast(multiplier);

  shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>> acaOp =
      DiscreteAcaBoundaryOperator<ValueType>::castToAca(op);
  const int symmetry = acaOp->symmetry();
  if (imagPart(multiplier) != 0. && symmetry & HERMITIAN &&
      !(symmetry & SYMMETRIC))
    throw std::runtime_error(
        "scaledAcaOperator(): multiplication of non-symmetric Hermitian"
        "matrices by scalars with non-zero imaginary part is not "
        "supported yet");
  shared_ptr<const AhmedBemBlcluster> scaledBlockCluster =
      acaOp->blockCluster();
  blcluster *nonConstScaledBlockCluster =
      const_cast<blcluster *>(            // AHMED is not const-correct
          static_cast<const blcluster *>( // safe -- upcast
              scaledBlockCluster.get()));
  boost::shared_array<AhmedMblock *> scaledBlocks =
      allocateAhmedMblockArray<ValueType>(scaledBlockCluster.get());
  copyH(nonConstScaledBlockCluster, acaOp->blocks().get(), scaledBlocks.get());

  const size_t blockCount = scaledBlockCluster->nleaves();
  for (size_t b = 0; b < blockCount; ++b) {
    AhmedMblock *block = scaledBlocks[b];
    typename AhmedTypeTraits<ValueType>::Type *data = block->getdata();
    if (block->isLrM()) // low-rank
      for (size_t i = 0; i < block->getn1() * block->rank(); ++i)
        data[i] *= ahmedMultiplier;
    else {
      // we don't support complex Hermitian blocks yet
      assert(!(block->isHeM() && boost::is_complex<ValueType>()));
      if (!block->isLtM()) // dense, but not lower-triangular
        for (size_t i = 0; i < block->nvals(); ++i)
          data[i] *= ahmedMultiplier;
      else { // lower-triangular
        // diagonal entries have special meaning and should not be rescaled
        size_t size = block->getn1();
        assert(block->getn2() == size);
        for (size_t c = 0; c < size; ++c)
          for (size_t r = 1; r < size; ++r)
            data[c * size + r] *= ahmedMultiplier;
      }
    }
  }

  int scaledSymmetry = symmetry;
  if (imagPart(multiplier) != 0.)
    scaledSymmetry &= ~HERMITIAN;
  shared_ptr<const DiscreteBoundaryOperator<ValueType>> result(
      new DiscreteAcaBoundaryOperator<ValueType>(
          acaOp->rowCount(), acaOp->columnCount(), acaOp->eps(),
          acaOp->maximumRank(), scaledSymmetry, scaledBlockCluster,
          scaledBlocks, acaOp->domainPermutation(), acaOp->rangePermutation(),
          acaOp->parallelizationOptions()));
  return result;
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>> scaledAcaOperator(
    const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op,
    const ValueType &multiplier) {
  return scaledAcaOperator(multiplier, op);
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>>
acaOperatorApproximateLuInverse(
    const shared_ptr<const DiscreteBoundaryOperator<ValueType>> &op,
    double delta) {
  shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>> acaOp =
      DiscreteAcaBoundaryOperator<ValueType>::castToAca(op);
  shared_ptr<const DiscreteBoundaryOperator<ValueType>> result(
      new AcaApproximateLuInverse<ValueType>(*acaOp, delta));
  return result;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteAcaBoundaryOperator);

#define INSTANTIATE_FREE_FUNCTIONS(RESULT)                                     \
  template shared_ptr<const DiscreteBoundaryOperator<RESULT>> acaOperatorSum(  \
      const shared_ptr<const DiscreteBoundaryOperator<RESULT>> &op1,           \
      const shared_ptr<const DiscreteBoundaryOperator<RESULT>> &op2,           \
      double eps, int maximumRank);                                            \
  template shared_ptr<const DiscreteBoundaryOperator<RESULT>>                  \
  scaledAcaOperator(                                                           \
      const RESULT &multiplier,                                                \
      const shared_ptr<const DiscreteBoundaryOperator<RESULT>> &op);           \
  template shared_ptr<const DiscreteBoundaryOperator<RESULT>>                  \
  scaledAcaOperator(                                                           \
      const shared_ptr<const DiscreteBoundaryOperator<RESULT>> &op,            \
      const RESULT &multiplier);                                               \
  template shared_ptr<const DiscreteBoundaryOperator<RESULT>>                  \
  acaOperatorApproximateLuInverse(                                             \
      const shared_ptr<const DiscreteBoundaryOperator<RESULT>> &op,            \
      double delta)

#if defined(ENABLE_SINGLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(float);
#endif

#if defined(ENABLE_SINGLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_BASIS_FUNCTIONS) ||                                \
     defined(ENABLE_COMPLEX_KERNELS))
INSTANTIATE_FREE_FUNCTIONS(std::complex<float>);
#endif

#if defined(ENABLE_DOUBLE_PRECISION)
INSTANTIATE_FREE_FUNCTIONS(double);
#endif

#if defined(ENABLE_DOUBLE_PRECISION) &&                                        \
    (defined(ENABLE_COMPLEX_BASIS_FUNCTIONS) ||                                \
     defined(ENABLE_COMPLEX_KERNELS))
INSTANTIATE_FREE_FUNCTIONS(std::complex<double>);
#endif

} // namespace Bempp

#endif // WITH_AHMED
