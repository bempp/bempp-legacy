// Copyright (C) 2013 by the BEM++ Authors
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

#ifdef WITH_AHMED

#include "block_coalescer.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include "ahmed_aux.hpp"

// This is a workaround of the problem of the abs() function being declared
// both in Epetra and in AHMED. It relies of the implementation detail (!) that
// in Epetra the declaration of abs is put between #ifndef __IBMCPP__ ...
// #endif. So it may well break in future versions of Trilinos. The ideal
// solution would be for AHMED to use namespaces.
#ifndef __IBMCPP__
#define __IBMCPP__
#include <Epetra_CrsMatrix.h>
#undef __IBMCPP__
#else
#include <Epetra_CrsMatrix.h>
#endif

namespace Bempp {

namespace {

void getLeafBlclusters(blcluster *bc, blcluster **leafClusters,
                       unsigned leafCount) {
  if (!bc)
    throw std::runtime_error("bc null");
  if (bc->isleaf()) {
    if (bc->getidx() >= leafCount)
      throw std::runtime_error("idx too large");
    assert(bc->getidx() < leafCount);
    leafClusters[bc->getidx()] = bc;
  } else {
    for (unsigned rs = 0; rs < bc->getnrs(); ++rs)
      for (unsigned cs = 0; cs < bc->getncs(); ++cs)
        getLeafBlclusters(bc->getson(rs, cs), leafClusters, leafCount);
  }
}

} // namespace

template <typename ValueType>
BlockCoalescer<ValueType>::BlockCoalescer(
    blcluster *blclustersRoot, blcluster *decomposedBlclustersRoot,
    const shared_ptr<const Epetra_CrsMatrix> &permutedTestGlobalToFlatLocalMap,
    const shared_ptr<const Epetra_CrsMatrix> &permutedTrialGlobalToFlatLocalMap,
    const boost::shared_array<AhmedMblock *> &blocks,
    const boost::shared_array<AhmedMblock *> &decomposedBlocks,
    const AcaOptions &acaOptions)
    : m_permutedTestGlobalToFlatLocalMap(permutedTestGlobalToFlatLocalMap),
      m_permutedTrialGlobalToFlatLocalMap(permutedTrialGlobalToFlatLocalMap),
      m_blocks(blocks), m_decomposedBlocks(decomposedBlocks),
      m_acaOptions(acaOptions) {
  unsigned blclusterCount = blclustersRoot->nleaves();
  if (decomposedBlclustersRoot->nleaves() != blclusterCount)
    throw std::runtime_error("BlockCoalescer::BlockCoalescer(): "
                             "cluster trees do not match");
  m_blclusters.reset(new blcluster *[blclusterCount]);
  m_decomposedBlclusters.reset(new blcluster *[blclusterCount]);
  getLeafBlclusters(blclustersRoot, m_blclusters.get(), blclusterCount);
  getLeafBlclusters(decomposedBlclustersRoot, m_decomposedBlclusters.get(),
                    blclusterCount);
}

template <typename ValueType>
void BlockCoalescer<ValueType>::coalesceBlock(unsigned index) {
  // early exit?
  if (!m_permutedTestGlobalToFlatLocalMap &&
      !m_permutedTrialGlobalToFlatLocalMap) {
    delete m_blocks[index];
    m_blocks[index] = 0;
    std::swap(m_decomposedBlocks[index], m_blocks[index]);
    return;
  }

  if (m_decomposedBlocks[index]->isGeM()) {
    if (m_decomposedBlocks[index]->isLtM() ||
        m_decomposedBlocks[index]->isUtM() ||
        m_decomposedBlocks[index]->isSyM() ||
        m_decomposedBlocks[index]->isHeM())
      throw std::runtime_error("BlockCoalescer::coalesceBlock(): "
                               "triangular and symmetric blocks are not "
                               "supported yet");
    coalesceDenseBlock(index);
  } else
    coalesceLowRankBlock(index);
}

template <typename ValueType>
void BlockCoalescer<ValueType>::coalesceDenseBlock(unsigned index) {
  // todo: handle the case of null m_permuted* pointers
  const blcluster *bc = m_blclusters[index];
  const blcluster *dbc = m_decomposedBlclusters[index];
  const ValueType *decomposedData = m_decomposedBlocks[index]->getdata();
  const unsigned rowCount = bc->getn1();
  const unsigned colCount = bc->getn2();
  const unsigned decomposedRowCount = dbc->getn1();
  const unsigned decomposedColCount = dbc->getn2();

  // this will be deallocated by freembls
  m_blocks[index] = new mblock<ValueType>(rowCount, colCount);
  m_blocks[index]->init0_GeM(rowCount, colCount);
  ValueType *data = m_blocks[index]->getdata();

  double ONE = 1.;
  int colIndex;
  int trialEntryCount, testEntryCount;
  double *trialValues, *testValues;
  int *trialIndices, *testIndices;

  colIndex = bc->getb2();
  for (unsigned j = 0; j < decomposedColCount; ++j, ++colIndex) {
    if (m_permutedTrialGlobalToFlatLocalMap)
      m_permutedTrialGlobalToFlatLocalMap->ExtractMyRowView(
          dbc->getb2() + j, trialEntryCount, trialValues, trialIndices);
    else {
      trialEntryCount = 1;
      trialValues = &ONE;
      trialIndices = &colIndex;
    }
    for (int trialDof = 0; trialDof < trialEntryCount; ++trialDof) {
      int n = trialIndices[trialDof] - bc->getb2();
      if (n >= colCount)
        throw std::runtime_error("BlockCoalescer::coalesceDenseBlock(): "
                                 "block is not independent");
      int rowIndex = bc->getb1();
      for (unsigned i = 0; i < decomposedRowCount; ++i, ++rowIndex) {
        if (m_permutedTestGlobalToFlatLocalMap)
          m_permutedTestGlobalToFlatLocalMap->ExtractMyRowView(
              dbc->getb1() + i, testEntryCount, testValues, testIndices);
        else {
          testEntryCount = 1;
          testValues = &ONE;
          testIndices = &rowIndex;
        }
        for (int testDof = 0; testDof < testEntryCount; ++testDof) {
          int m = testIndices[testDof] - bc->getb1();
          if (m >= rowCount)
            throw std::runtime_error("BlockCoalescer::coalesceDenseBlock(): "
                                     "block is not independent");
          data[m + n * rowCount] += testValues[testDof] *
                                    trialValues[trialDof] *
                                    decomposedData[i + j * decomposedRowCount];
        }
      }
    }
  }
  delete m_decomposedBlocks[index];
  m_decomposedBlocks[index] = 0;
}

template <typename ValueType>
void BlockCoalescer<ValueType>::coalesceLowRankBlock(unsigned index) {
  const blcluster *bc = m_blclusters[index];
  const blcluster *dbc = m_decomposedBlclusters[index];
  ValueType *decomposedData = m_decomposedBlocks[index]->getdata();
  const unsigned rowCount = bc->getn1();
  const unsigned colCount = bc->getn2();
  const unsigned decomposedRowCount = dbc->getn1();
  const unsigned decomposedColCount = dbc->getn2();
  const unsigned rank = m_decomposedBlocks[index]->rank();

  ValueType *U = 0, *V = 0;
  boost::scoped_array<ValueType> new_U, new_V;

  // this will be deallocated by freembls
  m_blocks[index] = new mblock<ValueType>(rowCount, colCount);

  int trialEntryCount, testEntryCount;
  double *trialValues, *testValues;
  int *trialIndices, *testIndices;

  if (m_permutedTestGlobalToFlatLocalMap) {
    new_U.reset(new ValueType[rank * rowCount]);
    U = new_U.get();
    std::fill(U, U + rank * rowCount, static_cast<ValueType>(0));
    for (unsigned i = 0; i < decomposedRowCount; ++i) {
      m_permutedTestGlobalToFlatLocalMap->ExtractMyRowView(
          dbc->getb1() + i, testEntryCount, testValues, testIndices);
      for (int testDof = 0; testDof < testEntryCount; ++testDof) {
        int m = testIndices[testDof] - bc->getb1();
        if (m >= rowCount)
          throw std::runtime_error("BlockCoalescer::coalesceLowRankBlock(): "
                                   "block is not independent");
        for (int k = 0; k < rank; ++k) {
          if (i + k * decomposedRowCount >= m_decomposedBlocks[index]->nvals())
            throw std::runtime_error("inv idx");
          U[m + k * rowCount] +=
              testValues[testDof] * decomposedData[i + k * decomposedRowCount];
        }
      }
    }
  } else
    U = decomposedData;

  if (m_permutedTrialGlobalToFlatLocalMap) {
    new_V.reset(new ValueType[rank * colCount]);
    V = new_V.get();
    std::fill(V, V + rank * colCount, static_cast<ValueType>(0));
    for (unsigned j = 0; j < decomposedColCount; ++j) {
      m_permutedTrialGlobalToFlatLocalMap->ExtractMyRowView(
          dbc->getb2() + j, trialEntryCount, trialValues, trialIndices);
      for (int trialDof = 0; trialDof < trialEntryCount; ++trialDof) {
        int n = trialIndices[trialDof] - bc->getb2();
        if (n >= colCount)
          throw std::runtime_error("BlockCoalescer::coalesceDenseBlock(): "
                                   "block is not independent");
        for (int k = 0; k < rank; ++k)
          V[n + k * colCount] +=
              trialValues[trialDof] *
              decomposedData
                  [rank * decomposedRowCount + j + k * decomposedColCount];
      }
    }
  } else
    V = decomposedData + rank * decomposedRowCount;

  m_blocks[index]
      ->cpyLrM_cmpr(rank, U, rowCount, V, colCount, m_acaOptions.eps, rank);
  delete m_decomposedBlocks[index];
  m_decomposedBlocks[index] = 0;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(BlockCoalescer);

} // namespace Bempp

#endif // WITH_AHMED
