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

#include "cluster_construction_helper.hpp"

#include "ahmed_aux.hpp"
#include "assembly_options.hpp"
#include "index_permutation.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../space/space.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/boost_make_shared_fwd.hpp"

namespace Bempp
{

template <typename BasisFunctionType>
void ClusterConstructionHelper<BasisFunctionType>::constructBemCluster(
        const Space<BasisFunctionType>& space,
        bool indexWithGlobalDofs,
        const AcaOptions& acaOptions,
        shared_ptr<AhmedBemCluster>& cluster,
        shared_ptr<IndexPermutation>& o2p,
        shared_ptr<IndexPermutation>& p2o)
{
#ifdef WITH_AHMED
    const size_t dofCount = indexWithGlobalDofs ?
                space.globalDofCount() : space.flatLocalDofCount();

    std::vector<unsigned int> o2pDofs(dofCount);
    std::vector<unsigned int> p2oDofs(dofCount);
    for (size_t i = 0; i < dofCount; ++i)
        o2pDofs[i] = i;
    for (size_t i = 0; i < dofCount; ++i)
        p2oDofs[i] = i;

    std::vector<Point3D<CoordinateType> > dofCenters;
    if (indexWithGlobalDofs)
        space.getGlobalDofPositions(dofCenters);
    else
        space.getFlatLocalDofPositions(dofCenters);

    // Use static_cast to convert from a pointer to Point3D to a pointer to its
    // descendant AhmedDofWrapper, which does not contain any new data members,
    // but just one additional method (the two structs should therefore be
    // binary compatible)
    AhmedDofType* ahmedDofCenters =
            static_cast<AhmedDofType*>(&dofCenters[0]);

    // NOTE: Ahmed uses names "op_perm" and "po_perm", which
    // correspond to BEM++'s "p2o" and "o2p", NOT the other way round.
    cluster = boost::make_shared<AhmedBemCluster>(
                ahmedDofCenters, &p2oDofs[0],
                0, dofCount, acaOptions.maximumBlockSize);
    cluster->createClusterTree(
        acaOptions.minimumBlockSize,
        &p2oDofs[0], &o2pDofs[0]);
    // cluster_pca stores a pointer to the first element of the array
    // dofCenters, but it is only used during construction of the
    // cluster tree. Now that the is done, we can deallocate
    // dofCenters. However, first we clear the references to this
    // array stored in the cluster tree, to detect any attempts to
    // access it more easily.
    cluster->clearDofPointers();
    o2p = boost::make_shared<IndexPermutation>(o2pDofs);
    p2o = boost::make_shared<IndexPermutation>(p2oDofs);
    // as far as I understand, cluster_pca doesn't store references to
    // p2oDofs or o2pDofs, so these arrays can now be safely deallocated.
#else // without Ahmed
    throw std::runtime_error("constructBemCluster(): "
                             "AHMED not available. Recompile BEM++ "
                             "with the symbol WITH_AHMED defined.");
#endif // WITH_AHMED
}

template <typename BasisFunctionType>
void ClusterConstructionHelper<BasisFunctionType>::constructBemCluster(
        const arma::Mat<CoordinateType>& points,
        const AcaOptions& acaOptions,
        shared_ptr<AhmedBemCluster>& cluster,
        shared_ptr<IndexPermutation>& o2p,
        shared_ptr<IndexPermutation>& p2o)
{
#ifdef WITH_AHMED
    const size_t pointCount = points.n_cols;
    const int dim = points.n_rows;
    if (dim > 3)
        throw std::invalid_argument(
                "ClusterConstructionHelper::constructBemCluster(): "
                "points from the array 'points' must have at most 3 coordinates");

    std::vector<unsigned int> o2pPoints(pointCount);
    std::vector<unsigned int> p2oPoints(pointCount);
    for (size_t i = 0; i < pointCount; ++i)
        o2pPoints[i] = i;
    for (size_t i = 0; i < pointCount; ++i)
        p2oPoints[i] = i;

    std::vector<Point3D<CoordinateType> > dofCenters;
    dofCenters.resize(pointCount);
    for (size_t p = 0; p < pointCount; ++p) {
        dofCenters[p].x = (dim >= 1) ? points(0, p) : 0.;
        dofCenters[p].y = (dim >= 2) ? points(1, p) : 0.;
        dofCenters[p].z = (dim >= 3) ? points(2, p) : 0.;
    }

    // Use static_cast to convert from a pointer to Point3D to a pointer to its
    // descendant AhmedDofWrapper, which does not contain any new data members,
    // but just one additional method (the two structs should therefore be
    // binary compatible)
    AhmedDofType* ahmedDofCenters =
            static_cast<AhmedDofType*>(&dofCenters[0]);

    // NOTE: Ahmed uses names "op_perm" and "po_perm", which
    // correspond to BEM++'s "p2o" and "o2p", NOT the other way round.
    cluster = boost::make_shared<AhmedBemCluster>(
                ahmedDofCenters, &p2oPoints[0],
                0, pointCount, acaOptions.maximumBlockSize);
    cluster->createClusterTree(
        acaOptions.minimumBlockSize,
        &p2oPoints[0], &o2pPoints[0]);
    // cluster_pca stores a pointer to the first element of the array
    // dofCenters, but it is only used during construction of the
    // cluster tree. Now that the is done, we can deallocate
    // dofCenters. However, first we clear the references to this
    // array stored in the cluster tree, to detect any attempts to
    // access it more easily.
    cluster->clearDofPointers();
    o2p = boost::make_shared<IndexPermutation>(o2pPoints);
    p2o = boost::make_shared<IndexPermutation>(p2oPoints);
    // as far as I understand, cluster_pca doesn't store references to
    // p2oDofs or o2pDofs, so these arrays can now be safely deallocated.
#else // without Ahmed
    throw std::runtime_error("constructBemCluster(): "
                             "AHMED not available. Recompile BEM++ "
                             "with the symbol WITH_AHMED defined.");
#endif // WITH_AHMED
}


template <typename BasisFunctionType>
std::auto_ptr<typename ClusterConstructionHelper<
                  BasisFunctionType>::AhmedBemBlcluster>
ClusterConstructionHelper<BasisFunctionType>::constructBemBlockCluster(
        const AcaOptions& acaOptions,
        bool symmetric,
        /* input parameter (effectively const */
        AhmedBemCluster& testCluster,
        /* input parameter (effectively const */
        AhmedBemCluster& trialCluster,
        /* output parameter */
        unsigned int& blockCount)
{
#ifdef WITH_AHMED
    std::auto_ptr<AhmedBemBlcluster> blockCluster(
                new AhmedBemBlcluster(&testCluster, &trialCluster));
    blockCount = 0;
    if (symmetric)
        blockCluster->subdivide_sym(&testCluster,
                                    acaOptions.eta * acaOptions.eta,
                                    blockCount);
    else
        blockCluster->subdivide(&testCluster, &trialCluster,
                                acaOptions.eta * acaOptions.eta,
                                blockCount);
    assert(blockCount == blockCluster->nleaves());
    return blockCluster;
#else // without Ahmed
    throw std::runtime_error("constructBlockBemCluster(): "
                             "AHMED not available. Recompile BEM++ "
                             "with the symbol WITH_AHMED defined.");
#endif // WITH_AHMED
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(ClusterConstructionHelper);

} // namespace Bempp
