#include "ahmed_leaf_cluster_array.hpp"

#ifdef WITH_AHMED

#include <algorithm>
#include <cmath>

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
    m_size(0) {
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

} // namespace Bempp

#endif
