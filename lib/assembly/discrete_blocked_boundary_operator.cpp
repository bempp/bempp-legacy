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

#include "bempp/common/config_trilinos.hpp"
#include "bempp/common/config_ahmed.hpp"
#include "../common/eigen_support.hpp"

#include "discrete_blocked_boundary_operator.hpp"

#include "discrete_aca_boundary_operator.hpp"
#ifdef WITH_AHMED
#include "ahmed_aux.hpp"
#endif
#include "../common/to_string.hpp"
#include "../fiber/_4d_array.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <numeric>
#ifdef WITH_TRILINOS
#include <Thyra_DefaultSpmdVectorSpace_decl.hpp>
#endif // WITH_TRILINOS

namespace Bempp {

// Functions used in
// DiscreteBlockBoundaryOperator::asDiscreteAcaBoundaryOperator().
namespace {

template <typename T> void dump(const Fiber::_2dArray<T> &a) {
  for (int i = 0; i < a.extent(0); ++i) {
    for (int j = 0; j < a.extent(1); ++j)
      std::cout << a(i, j) << "\t";
    std::cout << std::endl;
  }
}

template <typename T> void dump(const std::vector<T> &a) {
  for (int i = 0; i < a.size(); ++i)
    std::cout << a[i] << std::endl;
}

#ifdef WITH_AHMED
template <typename ValueType>
void collectMblocks(
    const Fiber::_2dArray<
        shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>>> &acaBlocks,
    const std::vector<size_t> &rowCounts,
    const std::vector<size_t> &columnCounts,
    std::vector<typename DiscreteAcaBoundaryOperator<
        ValueType>::AhmedConstMblockArray> &allSharedMblocks,
    typename DiscreteAcaBoundaryOperator<ValueType>::AhmedMblockArray &mblocks,
    Fiber::_2dArray<size_t> &indexOffsets, size_t &totalMblockCount) {
  typedef DiscreteAcaBoundaryOperator<ValueType> AcaOp;
  typedef typename AcaOp::AhmedMblock AhmedMblock;
  typedef typename AcaOp::AhmedMblockArray AhmedMblockArray;
  typedef typename AcaOp::AhmedConstMblockArray AhmedConstMblockArray;

  const size_t blockRowCount = acaBlocks.extent(0);
  const size_t blockColCount = acaBlocks.extent(1);

  // Put:
  // - into explSharedMblocks: arrays of mblocks making up the H-matrix
  //   of each block (including the single-element mblock arrays created for
  //   each empty block)
  // - into allSharedMblocks: the above plus the arrays of mblocks on which
  //   the non-empty blocks depend implicitly
  // - into mblockCounts: numbers of mblocks making up the H-matrix of each
  //   block
  // - into indexOffsets: numbers by which the indices of the mblocks making
  //   up the H-matrix of each block should be offset
  // - into totalMblockCount: the total size of all arrays from
  //   explSharedBlocks

  Fiber::_2dArray<AhmedConstMblockArray> explSharedMblocks(blockRowCount,
                                                           blockColCount);
  Fiber::_2dArray<size_t> mblockCounts(blockRowCount, blockColCount);
  indexOffsets.set_size(blockRowCount, blockColCount);
  std::fill(indexOffsets.begin(), indexOffsets.end(), 0);
  totalMblockCount = 0;
  for (size_t col = 0; col < blockColCount; ++col)
    for (size_t row = 0; row < blockRowCount; ++row) {
      if (acaBlocks(row, col)) {
        const AcaOp &activeOp = *acaBlocks(row, col);
        explSharedMblocks(row, col) = activeOp.blocks();
        allSharedMblocks.insert(allSharedMblocks.end(),
                                activeOp.sharedBlocks().begin(),
                                activeOp.sharedBlocks().end());
        mblockCounts(row, col) = activeOp.blockCount();
      } else {
        // Allocate a single empty block
        AhmedMblockArray dummyBlocks = allocateAhmedMblockArray<ValueType>(1);
        dummyBlocks[0] = new AhmedMblock(rowCounts[row], columnCounts[col]);
        explSharedMblocks(row, col) = dummyBlocks;
        mblockCounts(row, col) = 1;
      }
      allSharedMblocks.push_back(explSharedMblocks(row, col));
      indexOffsets(row, col) = totalMblockCount;
      totalMblockCount += mblockCounts(row, col);
    }

  // Create a shared array of pointers to all mblocks. Note: the array owns
  // only pointers, not mblocks themselves.
  mblocks.reset(new AhmedMblock *[totalMblockCount]);
  for (size_t col = 0; col < blockColCount; ++col)
    for (size_t row = 0; row < blockRowCount; ++row)
      std::copy(&explSharedMblocks(row, col)[0],
                &explSharedMblocks(row, col)[mblockCounts(row, col)],
                &mblocks[indexOffsets(row, col)]);
}

template <typename ValueType>
void copySonsAdjustingIndices(const blcluster *source, blcluster *dest,
                              // size_t rowOffset,
                              // size_t columnOffset,
                              size_t indexOffset) {
  assert(source);
  assert(dest);
  typedef typename DiscreteAcaBoundaryOperator<ValueType>::AhmedBemBlcluster
  AhmedBemBlcluster;

  unsigned int sonCount = source->getns();
  if (sonCount > 0) {
    // how I miss auto_array...
    std::vector<blcluster *> destSons(sonCount, 0);
    try {
      for (unsigned int row = 0, i = 0; row < source->getnrs(); ++row)
        for (unsigned int col = 0; col < source->getncs(); ++col, ++i) {
          const blcluster *sourceSon = source->getson(row, col);
          if (!sourceSon)
            continue;
          std::unique_ptr<blcluster> destSon;
          if (const AhmedBemBlcluster *bemSourceSon =
                  dynamic_cast<const AhmedBemBlcluster *>(sourceSon))
            destSon.reset(new AhmedBemBlcluster(
                // rowOffset + sourceSon->getb1(),
                // columnOffset + sourceSon->getb2(),
                dest->getb1() + sourceSon->getb1() - source->getb1(),
                dest->getb2() + sourceSon->getb2() - source->getb2(),
                sourceSon->getn1(), sourceSon->getn2()));
          else
            destSon.reset(new blcluster(
                // rowOffset + sourceSon->getb1(),
                // columnOffset + sourceSon->getb2(),
                dest->getb1() + sourceSon->getb1() - source->getb1(),
                dest->getb2() + sourceSon->getb2() - source->getb2(),
                sourceSon->getn1(), sourceSon->getn2()));
          copySonsAdjustingIndices<ValueType>(sourceSon, destSon.get(),
                                              // rowOffset, columnOffset,
                                              indexOffset);
          destSons[i] = destSon.release();
        }
      dest->setsons(source->getnrs(), source->getncs(), &destSons[0]);
    }
    catch (...) {
      for (unsigned int i = 0; i < sonCount; ++i)
        delete destSons[i];
      throw; // rethrow
    }
  } else {
    dest->setidx(indexOffset + source->getidx());
    dest->setadm(dest->isadm());
    dest->setsep(dest->issep());
  }
}

template <typename ValueType>
double overallEps(const Fiber::_2dArray<
    shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>>> &acaBlocks) {
  double result = std::numeric_limits<double>::max();
  for (size_t col = 0; col < acaBlocks.extent(1); ++col)
    for (size_t row = 0; row < acaBlocks.extent(0); ++row)
      if (acaBlocks(row, col))
        result = std::min(result, acaBlocks(row, col)->eps());
  return result;
}

template <typename ValueType>
int overallMaximumRank(const Fiber::_2dArray<
    shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>>> &acaBlocks) {
  int result = 0;
  for (size_t col = 0; col < acaBlocks.extent(1); ++col)
    for (size_t row = 0; row < acaBlocks.extent(0); ++row)
      if (acaBlocks(row, col))
        result = std::max(result, acaBlocks(row, col)->maximumRank());
  return result;
}

template <typename ValueType>
int overallSymmetry(const Fiber::_2dArray<
    shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>>> &acaBlocks) {
  for (size_t col = 0; col < acaBlocks.extent(1); ++col)
    for (size_t row = 0; row < acaBlocks.extent(0); ++row)
      if (acaBlocks(row, col) && acaBlocks(row, col)->symmetry() != 0)
        throw std::runtime_error("DiscreteBlockedBoundaryOperator::"
                                 "asDiscreteAcaBoundaryOperator(): "
                                 "Joining of symmetric ACA operators "
                                 "is not supported yet");
  return 0;
}

template <typename ValueType>
Fiber::ParallelizationOptions
overallParallelizationOptions(const Fiber::_2dArray<
    shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>>> &acaBlocks) {
  for (size_t col = 0; col < acaBlocks.extent(1); ++col)
    for (size_t row = 0; row < acaBlocks.extent(0); ++row)
      if (acaBlocks(row, col))
        return acaBlocks(row, col)->parallelizationOptions();
  return Fiber::ParallelizationOptions(); // actually should never get here
}

template <typename ValueType>
void getDomain_p2os(
    const Fiber::_2dArray<
        shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>>> &acaBlocks,
    std::vector<std::vector<unsigned>> &p2os, size_t &totalColCount) {
  const size_t rowBlockCount = acaBlocks.extent(0);
  const size_t colBlockCount = acaBlocks.extent(1);
  p2os.clear();
  p2os.reserve(colBlockCount);
  totalColCount = 0;
  for (size_t cb = 0; cb < colBlockCount; ++cb) {
    int chosenRow = -1;
    for (size_t rb = 0; rb < rowBlockCount; ++rb)
      if (acaBlocks(rb, cb)) {
        if (chosenRow < 0)
          chosenRow = rb;
        else { // ensure that all permutations are compatible
          if (acaBlocks(chosenRow, cb)->domainPermutation() !=
              acaBlocks(rb, cb)->domainPermutation())
            throw std::runtime_error("getDomain_p2os(): "
                                     "Domain index permutations of blocks (" +
                                     toString(chosenRow) + ", " + toString(cb) +
                                     ") and (" + toString(rb) + ", " +
                                     toString(cb) + ") are not equal");
        }
      }
    if (chosenRow < 0)
      throw std::runtime_error("getDomain_p2os(): "
                               "Column " +
                               toString(cb) + " contains empty blocks only");
    p2os.push_back(
        acaBlocks(chosenRow, cb)->domainPermutation().unpermutedIndices());
    const size_t localSize = p2os.back().size();
    for (size_t i = 0; i < localSize; ++i)
      p2os.back()[i] += totalColCount;
    totalColCount += localSize;
  }
}

template <typename ValueType>
void getRange_p2os(
    const Fiber::_2dArray<
        shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>>> &acaBlocks,
    std::vector<std::vector<unsigned>> &p2os, size_t &totalRowCount) {
  const size_t rowBlockCount = acaBlocks.extent(0);
  const size_t colBlockCount = acaBlocks.extent(1);
  p2os.clear();
  p2os.reserve(rowBlockCount);
  totalRowCount = 0;
  for (size_t rb = 0; rb < rowBlockCount; ++rb) {
    int chosenCol = -1;
    for (size_t cb = 0; cb < colBlockCount; ++cb)
      if (acaBlocks(rb, cb)) {
        if (chosenCol < 0)
          chosenCol = cb;
        else { // ensure that all permutations are compatible
          if (acaBlocks(rb, chosenCol)->rangePermutation() !=
              acaBlocks(rb, cb)->rangePermutation())
            throw std::runtime_error("getRange_p2os(): "
                                     "Range index permutations of blocks (" +
                                     toString(rb) + ", " + toString(chosenCol) +
                                     ") and (" + toString(rb) + ", " +
                                     toString(cb) + ") are not equal");
        }
      }
    if (chosenCol < 0)
      throw std::runtime_error("getRange_p2os(): "
                               "Row " +
                               toString(rb) + " contains empty blocks only");
    p2os.push_back(
        acaBlocks(rb, chosenCol)->rangePermutation().unpermutedIndices());
    const size_t localSize = p2os.back().size();
    for (size_t i = 0; i < localSize; ++i)
      p2os.back()[i] += totalRowCount;
    totalRowCount += localSize;
  }
}

template <typename ValueType>
void getUniformSonSizes(
    const Fiber::_2dArray<
        shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>>> &acaBlocks,
    // rowSonSizes[i](j, k): number of dofs in row son #k of block #j at level
    // #i
    std::vector<Fiber::_2dArray<unsigned>> &rowSonSizes,
    std::vector<Fiber::_2dArray<unsigned>> &colSonSizes) {
  const size_t rowBlockCount = acaBlocks.extent(0);
  const size_t colBlockCount = acaBlocks.extent(1);

  Fiber::_2dArray<const blcluster *> clusters(rowBlockCount, colBlockCount);
  bool anyBlocksAreNull = false;
  for (size_t col = 0; col < colBlockCount; ++col)
    for (size_t row = 0; row < rowBlockCount; ++row)
      if (acaBlocks(row, col))
        clusters(row, col) = acaBlocks(row, col)->blockCluster().get();
      else
        anyBlocksAreNull = true;

  if (anyBlocksAreNull) {
    // std::cout << "Null block" << std::endl;
    return;
  }

  Fiber::_4dArray<const blcluster *> parentClusters(rowBlockCount, 1,
                                                    colBlockCount, 1);
  for (size_t cb = 0; cb < colBlockCount; ++cb)
    for (size_t rb = 0; rb < rowBlockCount; ++rb)
      parentClusters(rb, 0, cb, 0) = clusters(rb, cb);
  size_t rowSonCount = 1;
  size_t colSonCount = 1;
  while (true) {
    size_t rowParentCount = parentClusters.extent(1);
    size_t colParentCount = parentClusters.extent(3);
    rowSonCount = parentClusters(0, 0, 0, 0)->getnrs();
    colSonCount = parentClusters(0, 0, 0, 0)->getncs();
    if (rowSonCount * colSonCount < 2) {
      // std::cout << "Leaf node" << std::endl;
      return;
    }

    for (size_t cb = 0; cb < colBlockCount; ++cb)
      for (size_t rb = 0; rb < rowBlockCount; ++rb)
        for (size_t cp = 0; cp < colParentCount; ++cp)
          for (size_t rp = 0; rp < rowParentCount; ++rp)
            if (parentClusters(rb, rp, cb, cp)->getnrs() != rowSonCount ||
                parentClusters(rb, rp, cb, cp)->getncs() != colSonCount) {
              // std::cout << "Mismatch in number of sons" << std::endl;
              return;
            }
    Fiber::_2dArray<unsigned> newRowSonSizes(rowBlockCount,
                                             rowParentCount * rowSonCount);
    Fiber::_2dArray<unsigned> newColSonSizes(colBlockCount,
                                             colParentCount * colSonCount);
    Fiber::_4dArray<const blcluster *> newClusters(
        rowBlockCount, rowSonCount * rowParentCount, colBlockCount,
        colSonCount * colParentCount);
    for (size_t cb = 0; cb < colBlockCount; ++cb)
      for (size_t rb = 0; rb < rowBlockCount; ++rb)
        for (size_t cp = 0; cp < colParentCount; ++cp)
          for (size_t rp = 0; rp < rowParentCount; ++rp)
            for (size_t cs = 0; cs < colSonCount; ++cs)
              for (size_t rs = 0; rs < rowSonCount; ++rs) {
                size_t R = rp * rowSonCount + rs;
                size_t C = cp * colSonCount + cs;
                const blcluster *son =
                    parentClusters(rb, rp, cb, cp)->getson(rs, cs);
                if (rp == 0 && rs == 0)
                  newColSonSizes(cb, C) = son->getn2();
                else if (newColSonSizes(cb, C) != son->getn2()) {
                  // std::cout << "Mismatch in column son sizes" << std::endl;
                  return;
                }
                if (cp == 0 && cs == 0)
                  newRowSonSizes(rb, R) = son->getn1();
                else if (newRowSonSizes(rb, R) != son->getn1()) {
                  // std::cout << "Mismatch in row son sizes" << std::endl;
                  return;
                }
                newClusters(rb, R, cb, C) = son;
              }
    rowSonSizes.push_back(newRowSonSizes);
    colSonSizes.push_back(newColSonSizes);
    parentClusters = newClusters;
  }
}

enum SpaceType {
  DOMAIN_SPACE,
  RANGE_SPACE
};

template <typename ValueType>
std::vector<unsigned int> overall_o2p(
    SpaceType mode,
    const Fiber::_2dArray<
        shared_ptr<const DiscreteAcaBoundaryOperator<ValueType>>> &acaBlocks,
    const Fiber::_2dArray<unsigned> &deepestSonSizes
    // (i, j) -> number of columns in deepest son #j of column block #i
    ) {
  assert(mode == DOMAIN_SPACE || mode == RANGE_SPACE);
  bool domainMode = (mode == DOMAIN_SPACE);
  const size_t blockCount = acaBlocks.extent(domainMode ? 1 : 0);
  assert(deepestSonSizes.extent(0) == blockCount);
  const size_t deepestSonCount = deepestSonSizes.extent(1);

  std::vector<std::vector<unsigned>> block_p2os;
  size_t totalDofCount;
  if (domainMode)
    getDomain_p2os(acaBlocks, block_p2os, totalDofCount /*, colBlockOffsets*/);
  else
    getRange_p2os(acaBlocks, block_p2os, totalDofCount /*, colBlockOffsets*/);
  assert(block_p2os.size() == blockCount);

  std::vector<unsigned> p2o;
  p2o.reserve(totalDofCount);

  std::vector<unsigned> bookmarks(blockCount, 0);
  for (unsigned s = 0; s < deepestSonCount; ++s)
    for (unsigned b = 0; b < blockCount; ++b)
      for (unsigned i = 0; i < deepestSonSizes(b, s); ++i)
        p2o.push_back(block_p2os[b][bookmarks[b]++]);

  std::vector<unsigned> o2p(totalDofCount);
  for (size_t i = 0; i < totalDofCount; ++i)
    o2p[p2o[i]] = i;
  return o2p;
}

void checkConsistency(const blcluster *b) {
  if (b->getnrs() == 0 || b->getncs() == 0)
    return;
  unsigned b1 = b->getb1();
  for (unsigned r = 0; r < b->getnrs(); ++r) {
    unsigned b2 = b->getb2();
    for (unsigned c = 0; c < b->getncs(); ++c) {
      assert(b->getson(r, c));
      assert(b->getson(r, c)->getb1() == b1);
      assert(b->getson(r, c)->getb2() == b2);
      checkConsistency(b->getson(r, c));
      b2 += b->getson(r, c)->getn2();
    }
    assert(b2 == b->getb2() + b->getn2());
    b1 += b->getson(r, 0)->getn1();
  }
  assert(b1 == b->getb1() + b->getn1());
}
#endif // WITH_AHMED

} // namespace

template <typename ValueType>
DiscreteBlockedBoundaryOperator<ValueType>::DiscreteBlockedBoundaryOperator(
    const Fiber::_2dArray<shared_ptr<const Base>> &blocks,
    const std::vector<size_t> &rowCounts,
    const std::vector<size_t> &columnCounts)
    : m_blocks(blocks), m_rowCounts(rowCounts), m_columnCounts(columnCounts) {
  if (rowCounts.size() != blocks.extent(0))
    throw std::invalid_argument(
        "DiscreteBlockedBoundaryOperator::DiscreteBlockedBoundaryOperator(): "
        "rowCounts has incorrect length");
  if (columnCounts.size() != blocks.extent(1))
    throw std::invalid_argument(
        "DiscreteBlockedBoundaryOperator::DiscreteBlockedBoundaryOperator(): "
        "columnCounts has incorrect length");

  for (size_t col = 0; col < blocks.extent(1); ++col)
    for (size_t row = 0; row < blocks.extent(0); ++row)
      if (blocks(row, col)) {
        if (blocks(row, col)->rowCount() != rowCounts[row] ||
            blocks(row, col)->columnCount() != columnCounts[col])
          throw std::invalid_argument(
              "DiscreteBlockedBoundaryOperator::"
              "DiscreteBlockedBoundaryOperator(): "
              "block (" +
              toString(row) + ", " + toString(col) +
              ") "
              "has incorrect dimensions (" +
              toString(blocks(row, col)->rowCount()) + ", " +
              toString(blocks(row, col)->columnCount()) + "); expected (" +
              toString(rowCounts[row]) + ", " + toString(columnCounts[col]) +
              ")");
      }
#ifdef WITH_TRILINOS
  m_domainSpace = Thyra::defaultSpmdVectorSpace<ValueType>(
      std::accumulate(m_columnCounts.begin(), m_columnCounts.end(), 0));
  m_rangeSpace = Thyra::defaultSpmdVectorSpace<ValueType>(
      std::accumulate(m_rowCounts.begin(), m_rowCounts.end(), 0));
#endif
}

template <typename ValueType>
unsigned int DiscreteBlockedBoundaryOperator<ValueType>::rowCount() const {
  return std::accumulate(m_rowCounts.begin(), m_rowCounts.end(), 0);
}

template <typename ValueType>
unsigned int DiscreteBlockedBoundaryOperator<ValueType>::columnCount() const {
  return std::accumulate(m_columnCounts.begin(), m_columnCounts.end(), 0);
}

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>>
DiscreteBlockedBoundaryOperator<ValueType>::getComponent(int row,
                                                         int col) const {
  return m_blocks(row, col);
}

template <typename ValueType>
void DiscreteBlockedBoundaryOperator<ValueType>::addBlock(
    const std::vector<int> &rows, const std::vector<int> &cols,
    const ValueType alpha, Matrix<ValueType> &block) const {
  throw std::runtime_error(
      "DiscreteBlockedBoundaryOperator::DiscreteBlockedBoundaryOperator(): "
      "addBlock: not implemented yet");
}

template <typename ValueType>
shared_ptr<const DiscreteBlockedBoundaryOperator<ValueType>>
DiscreteBlockedBoundaryOperator<
    ValueType>::asDiscreteAcaBlockedBoundaryOperator(double eps,
                                                     int maximumRank) const {
#ifdef WITH_AHMED
  size_t nrows = m_blocks.extent(0);
  size_t ncols = m_blocks.extent(1);
  Fiber::_2dArray<shared_ptr<const Base>> acaBlocks(nrows, ncols);
  for (size_t i = 0; i < nrows; i++) {
    for (size_t j = 0; j < ncols; j++) {
      if (m_blocks(i, j).get() != 0)
        acaBlocks(i, j) =
            m_blocks(i, j)->asDiscreteAcaBoundaryOperator(eps, maximumRank);
    }
  }
  return shared_ptr<const DiscreteBlockedBoundaryOperator<ValueType>>(
      new DiscreteBlockedBoundaryOperator(acaBlocks, m_rowCounts,
                                          m_columnCounts));
#else // WITH_AHMED
  throw std::runtime_error("DiscreteBlockedBoundaryOperator::"
                           "asDiscreteAcaBlockedBoundaryOperator(): "
                           "ACA operators are not supported because BEM++ "
                           "has been compiled without AHMED.");
#endif
}

#ifdef WITH_AHMED
// clusters is an array of block clusters lying in the same position in the
// trees
// of all the H matrices to be merged
template <typename ValueType>
void DiscreteBlockedBoundaryOperator<ValueType>::mergeHMatrices(
    unsigned currentLevel,
    const std::vector<Fiber::_2dArray<unsigned>> &rowSonSizes,
    const std::vector<Fiber::_2dArray<unsigned>> &colSonSizes,
    const Fiber::_2dArray<const blcluster *> clusters,
    const Fiber::_2dArray<size_t> indexOffsets, blcluster *result) const {
  typedef DiscreteAcaBoundaryOperator<ValueType> AcaOp;
  typedef typename AcaOp::AhmedBemBlcluster AhmedBemBlcluster;

  const size_t rowBlockCount = clusters.extent(0);
  const size_t colBlockCount = clusters.extent(1);
  assert(rowSonSizes.size() == colSonSizes.size());
  const size_t uniformLevelCount = rowSonSizes.size();
  if (currentLevel < uniformLevelCount) {
    // std::cout << "currentLevel = " << currentLevel << " < " <<
    // uniformLevelCount
    //           << std::endl;

    size_t rowSonCount = clusters(0, 0)->getnrs();
    size_t colSonCount = clusters(0, 0)->getncs();
    std::vector<blcluster *> branches(rowSonCount * colSonCount, 0);
    try {
      unsigned b1 = result->getb1();
      for (size_t rs = 0, i = 0; rs < rowSonCount; ++rs) {
        unsigned n1 = 0;
        for (size_t rb = 0; rb < rowBlockCount; ++rb)
          n1 += clusters(rb, 0)->getson(rs, 0)->getn1();
        unsigned b2 = result->getb2();
        for (size_t cs = 0; cs < colSonCount; ++cs, ++i) {
          unsigned n2 = 0;
          for (size_t cb = 0; cb < colBlockCount; ++cb)
            n2 += clusters(0, cb)->getson(0, cs)->getn2();
          std::unique_ptr<AhmedBemBlcluster> branch(
              new AhmedBemBlcluster(b1, b2, n1, n2));
          Fiber::_2dArray<const blcluster *> sons(rowBlockCount, colBlockCount);
          for (size_t rb = 0; rb < rowBlockCount; ++rb)
            for (size_t cb = 0; cb < colBlockCount; ++cb)
              sons(rb, cb) = clusters(rb, cb)->getson(rs, cs);

          mergeHMatrices(currentLevel + 1, rowSonSizes, colSonSizes, sons,
                         indexOffsets, branch.get());
          b2 += n2;
          branches[i] = branch.release();
        }
        b1 += n1;
      }
      result->setsons(rowSonCount, colSonCount, &branches[0]);
      checkConsistency(result);
    }
    catch (...) {
      for (size_t i = 0; i < branches.size(); ++i)
        delete branches[i];
      throw; // rethrow exception
    }
  } else if (currentLevel == uniformLevelCount) {
    // std::cout << "currentLevel = " << currentLevel << " == " <<
    // uniformLevelCount
    //           << std::endl;
    // Get number of elementary rows and columns in all rows and columns of
    // blocks
    std::vector<size_t> blockRowCount(rowBlockCount);
    std::vector<size_t> blockColCount(colBlockCount);

    if (uniformLevelCount == 0) {
      assert(rowBlockCount == m_rowCounts.size());
      assert(colBlockCount == m_columnCounts.size());
      blockRowCount = m_rowCounts;
      blockColCount = m_columnCounts;
    } else {
      for (size_t cb = 0; cb < colBlockCount; ++cb)
        for (size_t rb = 0; rb < rowBlockCount; ++rb) {
          const blcluster *curBlock = clusters(rb, cb);
          assert(curBlock);
          if (cb == 0)
            blockRowCount[rb] = curBlock->getn1();
          else
            assert(blockRowCount[rb] == curBlock->getn1());
          if (rb == 0)
            blockColCount[cb] = curBlock->getn2();
          else
            assert(blockColCount[cb] == curBlock->getn2());
        }
    }

    std::vector<blcluster *> branches(rowBlockCount * colBlockCount, 0);
    try {
      unsigned b1 = result->getb1();
      for (size_t rb = 0, i = 0; rb < rowBlockCount; ++rb) {
        unsigned b2 = result->getb2();
        for (size_t cb = 0; cb < colBlockCount; ++cb, ++i) {
          std::unique_ptr<AhmedBemBlcluster> branch(new AhmedBemBlcluster(
              b1, b2, blockRowCount[rb], blockColCount[cb]));
          if (clusters(rb, cb))
            copySonsAdjustingIndices<ValueType>(clusters(rb, cb), branch.get(),
                                                indexOffsets(rb, cb));
          else {
            branch->setidx(indexOffsets(rb, cb));
            branch->setadm(true); // unsure about it
            branch->setsep(true); // unsure about it
          }
          checkConsistency(branch.get());
          b2 += blockColCount[cb];
          branches[i] = branch.release();
        }
        b1 += blockRowCount[rb];
      }
      result->setsons(rowBlockCount, colBlockCount, &branches[0]);
      checkConsistency(result);
    }
    catch (...) {
      for (size_t i = 0; i < branches.size(); ++i)
        delete branches[i];
      throw; // rethrow exception
    }
  } else
    throw std::invalid_argument("mergeHMatrices(): invalid level");
}
#endif

template <typename ValueType>
shared_ptr<const DiscreteBoundaryOperator<ValueType>>
DiscreteBlockedBoundaryOperator<ValueType>::asDiscreteAcaBoundaryOperator(
    double eps, int maximumRank, bool interleave) const {
#ifdef WITH_AHMED
  // Get block counts
  const size_t rowBlockCount = m_blocks.extent(0);
  const size_t colBlockCount = m_blocks.extent(1);

  typedef DiscreteAcaBoundaryOperator<ValueType> AcaOp;
  typedef typename AcaOp::AhmedBemBlcluster AhmedBemBlcluster;
  typedef typename AcaOp::AhmedMblock AhmedMblock;
  typedef typename AcaOp::AhmedMblockArray AhmedMblockArray;
  typedef typename AcaOp::AhmedConstMblockArray AhmedConstMblockArray;

  // Convert blocks into ACA operators
  Fiber::_2dArray<shared_ptr<const AcaOp>> acaBlocks(rowBlockCount,
                                                     colBlockCount);
  for (size_t col = 0; col < colBlockCount; ++col)
    for (size_t row = 0; row < rowBlockCount; ++row)
      if (m_blocks(row, col))
        acaBlocks(row, col) = boost::shared_dynamic_cast<const AcaOp>(m_blocks(
            row, col)->asDiscreteAcaBoundaryOperator(eps, maximumRank));
  // else: perhaps create a new "empty" aca operator

  Fiber::_2dArray<const blcluster *> blockClusters(rowBlockCount,
                                                   colBlockCount);
  for (size_t col = 0; col < colBlockCount; ++col)
    for (size_t row = 0; row < rowBlockCount; ++row)
      if (acaBlocks(row, col))
        blockClusters(row, col) = acaBlocks(row, col)->blockCluster().get();
      else
        blockClusters(row, col) = 0;

  std::vector<Fiber::_2dArray<unsigned>> rowSonSizes, colSonSizes;
  if (interleave) {
    getUniformSonSizes(acaBlocks, rowSonSizes, colSonSizes);
    // for (size_t i = 0; i < rowSonSizes.size(); ++i) {
    //     std::cout << "Row son sizes, level " << i << ":\n";
    //     dump(rowSonSizes[i]);
    // }
    // for (size_t i = 0; i < colSonSizes.size(); ++i) {
    //     std::cout << "Col son sizes, level " << i << ":\n";
    //     dump(colSonSizes[i]);
    // }
  }

  std::vector<AhmedConstMblockArray> allSharedMblocks;
  AhmedMblockArray mblocks;
  Fiber::_2dArray<size_t> indexOffsets(rowBlockCount, colBlockCount);
  size_t totalMblockCount = 0;
  collectMblocks(acaBlocks, m_rowCounts, m_columnCounts, allSharedMblocks,
                 mblocks, indexOffsets, totalMblockCount);

  const size_t totalRowCount =
      std::accumulate(m_rowCounts.begin(), m_rowCounts.end(), 0);
  const size_t totalColCount =
      std::accumulate(m_columnCounts.begin(), m_columnCounts.end(), 0);

  // Recursively copy block cluster trees, adjusting indices
  shared_ptr<AhmedBemBlcluster> rootCluster(
      new AhmedBemBlcluster(0, 0, totalRowCount, totalColCount));
  mergeHMatrices(0 /* current level */, rowSonSizes, colSonSizes, blockClusters,
                 indexOffsets, rootCluster.get());

  Fiber::_2dArray<unsigned> deepestRowSonSizes;
  if (rowSonSizes.empty()) {
    deepestRowSonSizes.set_size(rowBlockCount, 1);
    for (size_t row = 0; row < rowBlockCount; ++row)
      deepestRowSonSizes(row, 0) = m_rowCounts[row];
  } else
    deepestRowSonSizes = rowSonSizes.back();
  Fiber::_2dArray<unsigned> deepestColSonSizes;
  if (colSonSizes.empty()) {
    deepestColSonSizes.set_size(colBlockCount, 1);
    for (size_t col = 0; col < colBlockCount; ++col)
      deepestColSonSizes(col, 0) = m_columnCounts[col];
  } else
    deepestColSonSizes = colSonSizes.back();
  // std::cout << "deepestRowSonSizes:\n"; dump(deepestRowSonSizes);
  // std::cout << "deepestColSonSizes:\n"; dump(deepestColSonSizes);
  std::vector<unsigned int> o2pDomain =
      overall_o2p(DOMAIN_SPACE, acaBlocks, deepestColSonSizes);
  std::vector<unsigned int> o2pRange =
      overall_o2p(RANGE_SPACE, acaBlocks, deepestRowSonSizes);

  // Gather remaining data necessary to create the combined ACA operator
  const double overallEps_ = overallEps(acaBlocks);
  const int overallMaximumRank_ = overallMaximumRank(acaBlocks);
  const int symmetry = overallSymmetry(acaBlocks);
  const Fiber::ParallelizationOptions parallelOptions =
      overallParallelizationOptions(acaBlocks);

  // std::ofstream os("aca-block.ps");
  // psoutputGeH(os, rootCluster.get(),
  //             std::max(totalRowCount, totalColCount), mblocks.get());
  // os.close();

  shared_ptr<const DiscreteBoundaryOperator<ValueType>> result(
      new AcaOp(totalRowCount, totalColCount, overallEps_, overallMaximumRank_,
                symmetry, rootCluster, mblocks, o2pDomain, o2pRange,
                parallelOptions, allSharedMblocks));
  return result;

#else // WITH_AHMED
  throw std::runtime_error("DiscreteBoundaryOperator::"
                           "asDiscreteAcaBoundaryOperator(): "
                           "ACA operators are not supported because BEM++ "
                           "has been compiled without AHMED.");
#endif
}

#ifdef WITH_TRILINOS
template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
DiscreteBlockedBoundaryOperator<ValueType>::domain() const {
  return m_domainSpace;
}

template <typename ValueType>
Teuchos::RCP<const Thyra::VectorSpaceBase<ValueType>>
DiscreteBlockedBoundaryOperator<ValueType>::range() const {
  return m_rangeSpace;
}

template <typename ValueType>
bool DiscreteBlockedBoundaryOperator<ValueType>::opSupportedImpl(
    Thyra::EOpTransp M_trans) const {
  for (size_t col = 0; col < m_blocks.extent(1); ++col)
    for (size_t row = 0; row < m_blocks.extent(0); ++row)
      if (m_blocks(row, col) && !m_blocks(row, col)->opSupported(M_trans))
        return false;
  return true;
}
#endif // WITH_TRILINOS

template <typename ValueType>
void DiscreteBlockedBoundaryOperator<ValueType>::applyBuiltInImpl(
    const TranspositionMode trans, const Vector<ValueType> &x_in,
    Vector<ValueType> &y_inout, const ValueType alpha,
    const ValueType beta) const {
  bool transpose = (trans == TRANSPOSE || trans == CONJUGATE_TRANSPOSE);
  size_t y_count = transpose ? m_columnCounts.size() : m_rowCounts.size();
  size_t x_count = transpose ? m_rowCounts.size() : m_columnCounts.size();

  for (int yi = 0, y_start = 0; yi < y_count; ++yi) {
    size_t y_chunk_size = transpose ? m_columnCounts[yi] : m_rowCounts[yi];
    Eigen::Map<Vector<ValueType>> y_chunk(y_inout.data()+y_start,y_chunk_size);
//    arma::Col<ValueType> y_chunk(&y_inout[y_start], y_chunk_size,
//                                 false /* copy_aux_mem */);
    for (int xi = 0, x_start = 0; xi < x_count; ++xi) {
      size_t x_chunk_size = transpose ? m_rowCounts[xi] : m_columnCounts[xi];
      shared_ptr<const Base> op =
          transpose ? m_blocks(xi, yi) : m_blocks(yi, xi);
      //            arma::Col<ValueType> x_chunk(&x_in[colStart],
      // m_columnCounts[col],
      //                                         false /* copy_aux_mem */);
      if (xi == 0) {
        // This branch ensures that the "y += beta * y" part is done
        if (op)
          //                    op->apply(trans, x_chunk, y_chunk, alpha, beta);
          op->apply(trans, x_in.segment(x_start, x_chunk_size),
                    y_chunk, alpha, beta);
        else {
          if (beta == static_cast<ValueType>(0.))
            y_chunk.setZero();
          else
            y_chunk *= beta;
        }
      } else if (op)
        //                    op->apply(trans, x_chunk, y_chunk, alpha, 1.);
        op->apply(trans, x_in.segment(x_start, x_chunk_size ),
                  y_chunk, alpha, 1.);
      x_start += x_chunk_size;
    }
    y_start += y_chunk_size;
  }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(DiscreteBlockedBoundaryOperator);

} // namespace Bempp
