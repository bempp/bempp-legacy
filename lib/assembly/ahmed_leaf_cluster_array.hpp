#ifndef bempp_ahmed_leaf_cluster_array_hpp
#define bempp_ahmed_leaf_cluster_array_hpp

#include "config_ahmed.hpp"

#ifdef WITH_AHMED

#include <boost/scoped_array.hpp>

class blcluster;

namespace Bempp
{

class AhmedLeafClusterArray
{
public:
    // The parameter should actually be const, but Ahmed's gen_BlSequence lacks
    // const-correctness
    explicit AhmedLeafClusterArray(blcluster* clusterTree);

    size_t size() const {
        return m_size;
    }

    blcluster* operator[] (size_t n) {
        return m_leafClusters[n];
    }

    const blcluster* operator[] (size_t n) const {
        return m_leafClusters[n];
    }

    /** \brief Sort cluster list, putting biggest clusters first. */
    void sortAccordingToClusterSize();

private:
    boost::scoped_array<blcluster*> m_leafClusters;
    size_t m_size;
};

} // namespace Bempp

#endif // WITH_AHMED

#endif
