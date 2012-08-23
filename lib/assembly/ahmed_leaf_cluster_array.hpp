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

#ifndef bempp_ahmed_leaf_cluster_array_hpp
#define bempp_ahmed_leaf_cluster_array_hpp

#include "../common/common.hpp"

#include "bempp/common/config_ahmed.hpp"

#ifdef WITH_AHMED

#include "../common/boost_scoped_array_fwd.hpp"

class blcluster;

namespace Bempp
{

/** \ingroup weak_form_assembly_internal
 *  \brief Encapsulation of an array of pointers to AHMED's blcluster objects.
 */
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
