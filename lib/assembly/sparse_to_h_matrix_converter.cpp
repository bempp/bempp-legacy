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
#if defined(WITH_TRILINOS) && defined(WITH_AHMED)

#include "sparse_to_h_matrix_converter.hpp"

#include "ahmed_aux.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace {

// Adapted from AHMED
void copyH(unsigned n, mblock<double> **A, mblock<dcomp> **B) {
  for (unsigned i = 0; i < n; ++i) {
    assert(A[i]);
    delete B[i];
    B[i] = new mblock<dcomp>(A[i]->getn1(), A[i]->getn2());
    assert(B[i]);

    if (A[i]->isLrM()) {
      const unsigned k = A[i]->rank();
      B[i]->setrank(k);
    } else { // mbl is dense
      assert(A[i]->isGeM());
      if (A[i]->isHeM())
        B[i]->setHeM();
      else if (A[i]->isLtM())
        B[i]->setLtM();
      else if (A[i]->isUtM())
        B[i]->setUtM();
      else
        B[i]->setGeM();
    }
    for (unsigned int j = 0; j < A[i]->nvals(); ++j)
      B[i]->getdata()[j] = A[i]->getdata()[j];
  }
}

// Adapted from AHMED
void copyH(unsigned n, mblock<double> **A, mblock<scomp> **B) {
  for (unsigned i = 0; i < n; ++i) {
    assert(A[i]);
    delete B[i];
    B[i] = new mblock<scomp>(A[i]->getn1(), A[i]->getn2());
    assert(B[i]);

    if (A[i]->isLrM()) {
      const unsigned k = A[i]->rank();
      B[i]->setrank(k);
    } else { // mbl is dense
      assert(A[i]->isGeM());
      if (A[i]->isHeM())
        B[i]->setHeM();
      else if (A[i]->isLtM())
        B[i]->setLtM();
      else if (A[i]->isUtM())
        B[i]->setUtM();
      else
        B[i]->setGeM();
    }
    for (unsigned int j = 0; j < A[i]->nvals(); ++j)
      B[i]->getdata()[j] = A[i]->getdata()[j];
  }
}

} // namespace

namespace Bempp {

// Convert an array of mblock<FromType> objects (in) to an array of
// mblock<ToType> objects (out).
template <typename FromType, typename ToType> struct MblockArrayConverter {
  static void convert(unsigned int n,
                      const boost::shared_array<mblock<FromType> *> &in,
                      boost::shared_array<mblock<ToType> *> &out) {
    out = allocateAhmedMblockArray<ToType>(n);
    copyH(n, in.get(), out.get());
  }
};

template <typename Type> struct MblockArrayConverter<Type, Type> {
  static void convert(unsigned int n,
                      const boost::shared_array<mblock<Type> *> &in,
                      boost::shared_array<mblock<Type> *> &out) {
    out = in;
  }
};

template <typename ValueType>
void SparseToHMatrixConverter<ValueType>::constructHMatrix(
    int *rowOffsets, int *colIndices, double *values,
    std::vector<unsigned int> &domain_o2p, std::vector<unsigned int> &range_p2o,
    double eps, bbxbemblcluster<AhmedDofType, AhmedDofType> *blockCluster,
    boost::shared_array<AhmedMblock *> &mblocks, int &maximumRank) {
  typedef mblock<AhmedTypeTraits<double>::Type> AhmedDoubleMblock;
  AhmedDoubleMblock **rawDoubleMblocks = 0;
  const unsigned int blockCount = blockCluster->nleaves();
  convCRS_toGeH(
      values, reinterpret_cast<unsigned int *>(colIndices), // int* -> uint*
      reinterpret_cast<unsigned int *>(rowOffsets),         // int* -> uint*
      &domain_o2p[0], &range_p2o[0], eps, blockCluster, rawDoubleMblocks);
  boost::shared_array<AhmedDoubleMblock *> doubleMblocks(
      (rawDoubleMblocks), (AhmedMblockArrayDeleter(blockCount)));
  MblockArrayConverter<
      double, typename AhmedTypeTraits<ValueType>::Type>::convert(blockCount,
                                                                  doubleMblocks,
                                                                  mblocks);
  maximumRank = Hmax_rank(blockCluster, mblocks.get());
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(SparseToHMatrixConverter);

} // namespace Bempp

#endif // WITH_AHMED && WITH_TRILINOS
