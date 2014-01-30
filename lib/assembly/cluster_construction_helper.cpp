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

#ifdef WITH_AHMED
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

// Unused function. Avoids a warning.
// namespace
// {
// 
// void dumpBlockCluster(blcluster& bc, int indent)
// {
//     for (int i = 0; i < indent; ++i)
//         std::cout << " ";
//     std::cout << bc.isadm() << "\n";
//     for (int r = 0; r < bc.getnrs(); ++r)
//         for (int c = 0; c < bc.getncs(); ++c)
//             dumpBlockCluster(*bc.getson(r, c), indent + 4);
// }
// 
// } // namespace

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

    std::vector<BoundingBox<CoordinateType> > dofCenters;
    if (indexWithGlobalDofs)
        space.getGlobalDofBoundingBoxes(dofCenters);
    else
        space.getFlatLocalDofBoundingBoxes(dofCenters);

    // Use static_cast to convert from a pointer to BoundingBox to a pointer to
    // its descendant AhmedDofWrapper, which does not contain any new data
    // members, but just a few additional methods (the two structs should
    // therefore be binary compatible)
    AhmedDofType* ahmedDofCenters =
            static_cast<AhmedDofType*>(&dofCenters[0]);

    // NOTE: Ahmed uses names "op_perm" and "po_perm", which
    // correspond to BEM++'s "p2o" and "o2p", NOT the other way round.
//    std::cout << "making a cluster\n";
    bool strongAdmissibility = !indexWithGlobalDofs;
    cluster = boost::make_shared<AhmedBemCluster>(
                ahmedDofCenters, &p2oDofs[0],
                0, dofCount, acaOptions.maximumBlockSize,
                strongAdmissibility);
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
        int componentCount,
        const AcaOptions& acaOptions,
        shared_ptr<AhmedBemCluster>& cluster,
        shared_ptr<IndexPermutation>& o2p,
        shared_ptr<IndexPermutation>& p2o)
{
#ifdef WITH_AHMED
    const size_t pointCount = points.n_cols;
    const size_t dofCount = pointCount * componentCount;

    std::vector<unsigned int> o2pDofs(dofCount);
    std::vector<unsigned int> p2oDofs(dofCount);
    for (size_t i = 0; i < dofCount; ++i)
        o2pDofs[i] = i;
    for (size_t i = 0; i < dofCount; ++i)
        p2oDofs[i] = i;

    std::vector<BoundingBox<CoordinateType> > dofCenters;
    getComponentBoundingBoxes(points, componentCount, dofCenters);

    // Use static_cast to convert from a pointer to BoundingBox to a pointer to
    // its descendant AhmedDofWrapper, which does not contain any new data
    // members, but just a few additional methods (the two structs should
    // therefore be binary compatible)
    AhmedDofType* ahmedDofCenters =
            static_cast<AhmedDofType*>(&dofCenters[0]);

    // NOTE: Ahmed uses names "op_perm" and "po_perm", which
    // correspond to BEM++'s "p2o" and "o2p", NOT the other way round.
    cluster = boost::make_shared<AhmedBemCluster>(
                ahmedDofCenters, &p2oDofs[0],
                0, dofCount, acaOptions.maximumBlockSize);
    cluster->createClusterTree(
        acaOptions.minimumBlockSize, &p2oDofs[0], &o2pDofs[0]);

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
std::auto_ptr<typename ClusterConstructionHelper<
                  BasisFunctionType>::AhmedBemBlcluster>
ClusterConstructionHelper<BasisFunctionType>::constructBemBlockCluster(
        const AcaOptions& acaOptions,
        bool symmetric,
        AhmedBemCluster& testCluster,
        AhmedBemCluster& trialCluster,
        bool useStrongAdmissibilityCondition,
        unsigned int& blockCount)
{
#ifdef WITH_AHMED
    // Save original admissibility conditions used in clusters
    const bool testStrongAdm = testCluster.usingStrongAdmissibilityCondition();
    const bool trialStrongAdm = trialCluster.usingStrongAdmissibilityCondition();
    try {
        if (testStrongAdm != useStrongAdmissibilityCondition)
            testCluster.useStrongAdmissibilityCondition(
                        useStrongAdmissibilityCondition);
        if (trialStrongAdm != useStrongAdmissibilityCondition)
            trialCluster.useStrongAdmissibilityCondition(
                        useStrongAdmissibilityCondition);

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
        // Restore original admissibility condition
        if (testStrongAdm != useStrongAdmissibilityCondition)
            testCluster.useStrongAdmissibilityCondition(testStrongAdm);
        if (trialStrongAdm != useStrongAdmissibilityCondition)
            trialCluster.useStrongAdmissibilityCondition(trialStrongAdm);
        return blockCluster;
    } catch (...) {
        // Restore original admissibility condition
        if (testStrongAdm != useStrongAdmissibilityCondition)
            testCluster.useStrongAdmissibilityCondition(testStrongAdm);
        if (trialStrongAdm != useStrongAdmissibilityCondition)
            trialCluster.useStrongAdmissibilityCondition(trialStrongAdm);
        throw;
    }

#else // without Ahmed
    throw std::runtime_error("constructBlockBemCluster(): "
                             "AHMED not available. Recompile BEM++ "
                             "with the symbol WITH_AHMED defined.");
#endif // WITH_AHMED
}

template <typename BasisFunctionType>
void
ClusterConstructionHelper<BasisFunctionType>::truncateBemBlockCluster(
        blcluster* cluster, const blcluster* refCluster)
{
    if (refCluster->isleaf()) {
        if (!cluster->isleaf())
            cluster->setsons(0, 0, 0);
        cluster->setidx(refCluster->getidx());
    } else {
        if (cluster->getnrs() != refCluster->getnrs() ||
                cluster->getncs() != refCluster->getncs())
            throw std::runtime_error("truncateBemBlockCluster(): the cluster "
                                     "to be truncated is not a superset of the "
                                     "reference cluster");
        for (unsigned r = 0; r < refCluster->getnrs(); ++r)
            for (unsigned c = 0; c < refCluster->getncs(); ++c)
                truncateBemBlockCluster(cluster->getson(r, c),
                                        refCluster->getson(r, c));
    }
}

template <typename BasisFunctionType>
void
ClusterConstructionHelper<BasisFunctionType>::getComponentDofPositions(
        const arma::Mat<CoordinateType>& points,
        int componentCount,
        std::vector<Point3D<CoordinateType> >& positions)
{
    const size_t pointCount = points.n_cols;
    const int dim = points.n_rows;
    if (dim > 3)
        throw std::invalid_argument(
                "ClusterConstructionHelper::getComponentDofPositions(): "
                "points from the array 'points' must have at most 3 coordinates");

    positions.resize(pointCount * componentCount);
    for (size_t p = 0; p < pointCount; ++p)
        for (size_t c = 0; c < componentCount; ++c) {
            positions[p * componentCount + c].x = (dim >= 1) ? points(0, p) : 0.;
            positions[p * componentCount + c].y = (dim >= 2) ? points(1, p) : 0.;
            positions[p * componentCount + c].z = (dim >= 3) ? points(2, p) : 0.;
        }
}

template <typename BasisFunctionType>
void
ClusterConstructionHelper<BasisFunctionType>::getComponentBoundingBoxes(
        const arma::Mat<CoordinateType>& points,
        int componentCount,
        std::vector<BoundingBox<CoordinateType> >& boundingBoxes)
{
    const size_t pointCount = points.n_cols;
    const int dim = points.n_rows;
    if (dim > 3)
        throw std::invalid_argument(
                "ClusterConstructionHelper::getComponentDofBoundingBoxes(): "
                "points from the array 'points' must have at most 3 coordinates");

    boundingBoxes.resize(pointCount * componentCount);
    for (size_t p = 0; p < pointCount; ++p)
        for (size_t c = 0; c < componentCount; ++c) {
            BoundingBox<CoordinateType>& bb = boundingBoxes[p * componentCount + c];
            bb.reference.x = bb.lbound.x = bb.ubound.x = (dim >= 1) ? points(0, p) : 0.;
            bb.reference.y = bb.lbound.y = bb.ubound.y = (dim >= 2) ? points(1, p) : 0.;
            bb.reference.z = bb.lbound.z = bb.ubound.z = (dim >= 3) ? points(2, p) : 0.;
        }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS(ClusterConstructionHelper);

} // namespace Bempp

#endif // WITH_AHMED
