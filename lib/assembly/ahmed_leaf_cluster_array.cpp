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

#include "ahmed_leaf_cluster_array.hpp"

#ifdef WITH_AHMED

#include <algorithm>
#include <cmath>
#include <stdexcept>

#define BASMOD // prevent inclusion of Ahmed's basmod.h, which contains
               // a conflicting definition of swap()
using std::min;
using std::max;
inline unsigned pow2(unsigned k)
{
  return (1 << k);
}
#include <blcluster.h>
#include <bllist.h>

namespace Bempp
{

namespace
{

bool isFirstClusterBigger(const blcluster* cluster1, const blcluster* cluster2)
{
    return cluster1->getn1() * cluster1->getn2() >
            cluster2->getn1() * cluster2->getn2();
}

} // namespace

AhmedLeafClusterArray::AhmedLeafClusterArray(blcluster* clusterTree) :
    m_size(0)
{
    blcluster** leafClusters = 0;
    try {
        gen_BlSequence(clusterTree, leafClusters);
    }
    catch (...) {
        delete[] leafClusters;
        throw; // rethrow the exception
    }
    m_leafClusters.reset(leafClusters);
    m_size = clusterTree->nleaves();
}


void AhmedLeafClusterArray::sortAccordingToClusterSize()
{
    std::sort(&m_leafClusters[0], &m_leafClusters[m_size],
              isFirstClusterBigger);
}

void AhmedLeafClusterArray::startWithClusterOfIndex(size_t index)
{
    if (index >= m_size)
        throw std::invalid_argument("AhmedLeafClusterArray::"
                                    "startWithClusterOfIndex(): invalid index");
    for (size_t i = 0; i < m_size; ++i)
        if (m_leafClusters[i]->getidx() == index) {
            std::swap(m_leafClusters[0], m_leafClusters[i]);
            break;
        }
}

} // namespace Bempp

#endif
