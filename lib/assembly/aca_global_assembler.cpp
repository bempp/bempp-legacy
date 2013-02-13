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

#include "aca_global_assembler.hpp"

#include "assembly_options.hpp"
#include "cluster_construction_helper.hpp"
#include "evaluation_options.hpp"
#include "index_permutation.hpp"
#include "discrete_boundary_operator_composition.hpp"
#include "discrete_sparse_boundary_operator.hpp"

#include "../common/armadillo_fwd.hpp"
#include "../common/auto_timer.hpp"
#include "../common/boost_shared_array_fwd.hpp"
#include "../common/chunk_statistics.hpp"
#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_operators.hpp"
#include "../fiber/local_assembler_for_potential_operators.hpp"
#include "../fiber/serial_blas_region.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../space/space.hpp"

#include <stdexcept>
#include <fstream>
#include <iostream>

#include <boost/type_traits/is_complex.hpp>

#include <tbb/atomic.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_queue.h>

#ifdef WITH_AHMED
#include "ahmed_aux.hpp"

#ifdef __INTEL_COMPILER
#pragma warning(disable:381)
#endif

#include <apprx.h>

#ifdef __INTEL_COMPILER
#pragma warning(default:381)
#endif

#include "discrete_aca_boundary_operator.hpp"
#include "modified_aca.hpp"
#include "potential_operator_aca_assembly_helper.hpp"
#include "scattered_range.hpp"
#include "weak_form_aca_assembly_helper.hpp"
#endif

// #define DUMP_DENSE_BLOCKS // if defined, contents and DOF lists of blocks
                             // stored as dense matrices will be printed to the
                             // screen

namespace Bempp
{

// Body of parallel loop
namespace
{

#ifdef WITH_AHMED
template <typename BasisFunctionType, typename ResultType,
          typename AcaAssemblyHelper>
class AcaAssemblerLoopBody
{
    typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;
    typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
    typedef bemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;
    typedef mblock<typename AhmedTypeTraits<ResultType>::Type> AhmedMblock;
public:
    typedef tbb::concurrent_queue<size_t> LeafClusterIndexQueue;

    AcaAssemblerLoopBody(
            AcaAssemblyHelper& helper,
            AhmedLeafClusterArray& leafClusters,
            boost::shared_array<AhmedMblock*> blocks,
            const AcaOptions& options,
            tbb::atomic<size_t>& done,
            bool verbose,
            LeafClusterIndexQueue& leafClusterIndexQueue,
            bool symmetric,
            std::vector<ChunkStatistics>& stats) :
        m_helper(helper),
        m_leafClusters(leafClusters), m_blocks(blocks),
        m_options(options), m_done(done), m_verbose(verbose),
        m_leafClusterIndexQueue(leafClusterIndexQueue),
        m_symmetric(symmetric),
        m_stats(stats)
    {
    }

    template <typename Range>
    void operator() (const Range& r) const {
        const char* TEXT = "Approximating ... ";
        for (typename Range::const_iterator i = r.begin(); i != r.end(); ++i) {
            size_t leafClusterIndex = 0;
            if (!m_leafClusterIndexQueue.try_pop(leafClusterIndex)) {
                std::cerr << "AcaWeakFormAssemblerLoopBody::operator(): "
                             "Warning: try_pop failed; this shouldn't happen!"
                          << std::endl;
                continue;
            }
            m_stats[leafClusterIndex].valid = true;
            m_stats[leafClusterIndex].chunkStart = r.begin();
            m_stats[leafClusterIndex].chunkSize = r.size();
            m_stats[leafClusterIndex].startTime = tbb::tick_count::now();

            AhmedBemBlcluster* cluster =
                    dynamic_cast<AhmedBemBlcluster*>(m_leafClusters[leafClusterIndex]);
            if (m_symmetric)
                apprx_sym(m_helper, m_blocks[cluster->getidx()],
                          cluster, m_options.eps, m_options.maximumRank,
                          true /* complex_sym */);
            else {
                if (m_options.useAhmedAca)
                    apprx_unsym(m_helper, m_blocks[cluster->getidx()],
                                cluster, m_options.eps, m_options.maximumRank);
                else
                    apprx_unsym_shooting(
                                m_helper, m_blocks[cluster->getidx()],
                                cluster, m_options.eps, m_options.maximumRank);
            }
            m_stats[leafClusterIndex].endTime = tbb::tick_count::now();
            // TODO: recompress
            const int HASH_COUNT = 20;
            if (m_verbose)
                progressbar(std::cout, TEXT, (++m_done) - 1,
                            m_leafClusters.size(), HASH_COUNT, true);
        }

    }

private:
    AcaAssemblyHelper& m_helper;
    AhmedLeafClusterArray& m_leafClusters;
    boost::shared_array<AhmedMblock*> m_blocks;
    const AcaOptions& m_options;
    tbb::atomic<size_t>& m_done;
    bool m_verbose;
    LeafClusterIndexQueue& m_leafClusterIndexQueue;
    bool m_symmetric;
    std::vector<ChunkStatistics>& m_stats;
};

void reallyGetClusterIds(const cluster& clusterTree,
                         const std::vector<unsigned int>& p2oDofs,
                         std::vector<unsigned int>& clusterIds,
                         unsigned int& id)
{
    if (clusterTree.isleaf())
        for (unsigned int nDof = clusterTree.getnbeg(); nDof < clusterTree.getnend(); ++nDof)
            clusterIds[p2oDofs[nDof]] = id;
    else
        for (unsigned int nSon = 0; nSon < clusterTree.getns(); ++nSon)
            reallyGetClusterIds(*clusterTree.getson(nSon), p2oDofs, clusterIds, ++id);
}

void getClusterIds(const cluster& clusterTree,
                   const std::vector<unsigned int>& p2oDofs,
                   std::vector<unsigned int>& clusterIds)
{
    clusterIds.resize(p2oDofs.size());
    unsigned int id = 0;
    reallyGetClusterIds(clusterTree, p2oDofs, clusterIds, id);
}

template <typename T>
void save_arma_matrix(const arma::Mat<T>& a, const std::string& fname)
{
    std::ofstream out;
    out.precision(17);
    out.open(fname.c_str());
    arma::diskio::save_raw_ascii(a, out);
    out.close();
}
template <typename T>
void save_arma_matrix(const arma::Mat<std::complex<T> >& a, const std::string& fname)
{
    std::ofstream out;
    out.precision(17);
    out.open((fname + "-real.txt").c_str());
    arma::Mat<T> component = arma::real(a);
    arma::diskio::save_raw_ascii(component, out);
    out.close();
    out.open((fname + "-imag.txt").c_str());
    component = arma::imag(a);
    arma::diskio::save_raw_ascii(component, out);
    out.close();
}

template <typename ValueType>
void dumpDenseBlocks(
        typename DiscreteAcaBoundaryOperator<ValueType>::AhmedBemBlcluster* clusterTree,
        typename DiscreteAcaBoundaryOperator<ValueType>::AhmedMblockArray& blocks,
        const std::vector<unsigned int>& p2oRows,
        const std::vector<unsigned int>& p2oCols,
        const std::vector<Point3D<typename Fiber::ScalarTraits<ValueType>::RealType> >& rowDofs,
        const std::vector<Point3D<typename Fiber::ScalarTraits<ValueType>::RealType> >& colDofs)
{
    if (!clusterTree)
        return;
    typedef typename Fiber::ScalarTraits<ValueType>::RealType CoordinateType;
    typedef DiscreteAcaBoundaryOperator<ValueType> AcaOp;
    typedef typename AcaOp::AhmedDofType AhmedDofType;
    typedef typename AcaOp::AhmedMblock AhmedMblock;
    typedef bemcluster<AhmedDofType> Cluster;
    if (clusterTree->isleaf()) {
        unsigned int idx = clusterTree->getidx();
        std::cout << "LEAF; idx = " << idx << std::endl;
        if (clusterTree->isGeM(blocks.get()) && clusterTree->isadm()) {
            std::cout << "Dense block; " << clusterTree->getb1()
                      << " " << clusterTree->getb2()
                      << " " << clusterTree->getn1()
                      << " " << clusterTree->getn2() << "\n";
            // if (clusterTree->getn1() < 500 || clusterTree->getn2() < 500)
            //     return;
            Cluster* clRow = clusterTree->getcl1();
            assert(clRow);
            std::cout << "Row center of mass: ("
                      << clRow->getcom(0) << ", "
                      << clRow->getcom(1) << ", "
                      << clRow->getcom(2) << ")\n";
            std::cout << "Row icm: "
                      << clRow->geticom() - clusterTree->getb1() << std::endl;
            for (unsigned int nDof = clRow->getnbeg();
                 nDof < clRow->getnend(); ++nDof) {
                assert(nDof < p2oRows.size());
                assert(p2oRows[nDof] < rowDofs.size());
                const Point3D<CoordinateType> dofPos = rowDofs[p2oRows[nDof]];
                std::cout << "  Row dof #" << p2oRows[nDof] << " at "
                          << dofPos.x << ", " << dofPos.y << ", "
                          << dofPos.z << "\n";
            }
            Cluster* clCol = clusterTree->getcl2();
            assert(clCol);
            std::cout << "Col center of mass: ("
                      << clCol->getcom(0) << ", "
                      << clCol->getcom(1) << ", "
                      << clCol->getcom(2) << ")\n";
            std::cout << "Column icm: "
                      << clCol->geticom() - clusterTree->getb2() << std::endl;
            for (unsigned int nDof = clCol->getnbeg();
                 nDof < clCol->getnend(); ++nDof) {
                assert(nDof < p2oCols.size());
                assert(p2oCols[nDof] < colDofs.size());
                const Point3D<CoordinateType> dofPos = colDofs[p2oCols[nDof]];
                std::cout << "  Col dof #" << p2oCols[nDof] << " at "
                          << dofPos.x << ", " << dofPos.y << ", "
                          << dofPos.z << "\n";
            }
            AhmedMblock* block = blocks[idx];
            arma::Mat<ValueType> ablock(clusterTree->getn1(),
                                        clusterTree->getn2());
            for (size_t i = 0; i < block->nvals(); ++i)
                ablock[i] = block->getdata()[i];
            save_arma_matrix(ablock, "block-" + toString(idx) + ".txt");
        }
    }
    else
        for (unsigned int nRowSon = 0; nRowSon < clusterTree->getnrs(); ++nRowSon)
            for (unsigned int nColSon = 0; nColSon < clusterTree->getncs(); ++nColSon)
            dumpDenseBlocks<ValueType>(
                        dynamic_cast<typename AcaOp::AhmedBemBlcluster*>(
                            clusterTree->getson(nRowSon, nColSon)),
                        blocks,
                        p2oRows, p2oCols, rowDofs, colDofs);
}

template <typename AcaAssemblyHelper,
          typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteAcaBoundaryOperator<ResultType> >
assembleAcaOperator(
        AcaAssemblyHelper& helper,
        const shared_ptr<typename DiscreteAcaBoundaryOperator<ResultType>::
            AhmedBemBlcluster>& bemBlclusterTree,
        const ParallelizationOptions& parallelOptions,
        const AcaOptions& acaOptions,
        bool verbosityAtLeastDefault,
        bool symmetric,
        const shared_ptr<IndexPermutation>& test_o2pPermutation,
        const shared_ptr<IndexPermutation>& trial_o2pPermutation
#ifdef DUMP_DENSE_BLOCKS
        ,
        const shared_ptr<IndexPermutation>& test_p2oPermutation,
        const shared_ptr<IndexPermutation>& trial_p2oPermutation,
        const std::vector<Point3D<
            typename Fiber::ScalarTraits<ResultType>::RealType> >& testDofCenters,
        const std::vector<Point3D<
            typename Fiber::ScalarTraits<ResultType>::RealType> >& trialDofCenters
#endif // DUMP_DENSE_BLOCKS
        )
{
    typedef mblock<typename AhmedTypeTraits<ResultType>::Type> AhmedMblock;
    boost::shared_array<AhmedMblock*> blocks =
            allocateAhmedMblockArray<ResultType>(bemBlclusterTree.get());

    const size_t testDofCount = test_o2pPermutation->size();
    const size_t trialDofCount = trial_o2pPermutation->size();

    AhmedLeafClusterArray leafClusters(bemBlclusterTree.get());
    leafClusters.sortAccordingToClusterSize();
    const size_t leafClusterCount = leafClusters.size();

    int maxThreadCount = 1;
    if (!parallelOptions.isOpenClEnabled())
    {
        if (parallelOptions.maxThreadCount() == ParallelizationOptions::AUTO)
            maxThreadCount = tbb::task_scheduler_init::automatic;
        else
            maxThreadCount = parallelOptions.maxThreadCount();
    }
    tbb::task_scheduler_init scheduler(maxThreadCount);
    tbb::atomic<size_t> done;
    done = 0;

    std::vector<ChunkStatistics> chunkStats(leafClusterCount);

    typedef AcaAssemblerLoopBody<
            BasisFunctionType, ResultType, AcaAssemblyHelper> Body;
    typename Body::LeafClusterIndexQueue leafClusterIndexQueue;
    for (size_t i = 0; i < leafClusterCount; ++i)
        leafClusterIndexQueue.push(i);

    if (verbosityAtLeastDefault)
        std::cout << "About to start the ACA assembly loop" << std::endl;
    tbb::tick_count loopStart = tbb::tick_count::now();
    {
        Fiber::SerialBlasRegion region; // if possible, ensure that BLAS is single-threaded
        tbb::parallel_for(tbb::blocked_range<size_t>(0, leafClusterCount),
                          Body(helper, leafClusters, blocks, acaOptions, done,
                               verbosityAtLeastDefault,
                               leafClusterIndexQueue, symmetric, chunkStats));
    }
    tbb::tick_count loopEnd = tbb::tick_count::now();
    if (verbosityAtLeastDefault) {
        std::cout << "\n"; // the progress bar doesn't print the final \n
        std::cout << "ACA loop took " << (loopEnd - loopStart).seconds() << " s"
                  << std::endl;
    }

    // TODO: parallelise!
    if (acaOptions.recompress) {
        if (verbosityAtLeastDefault)
            std::cout << "About to start ACA agglomeration" << std::endl;
        agglH(bemBlclusterTree.get(), blocks.get(),
              acaOptions.eps, acaOptions.maximumRank);
        if (verbosityAtLeastDefault)
            std::cout << "Agglomeration finished" << std::endl;
    }

#ifdef DUMP_DENSE_BLOCKS
    dumpDenseBlocks<ResultType>(bemBlclusterTree.get(), blocks,
                                test_p2oPermutation->permutedIndices(),
                                trial_p2oPermutation->permutedIndices(),
                                testDofCenters, trialDofCenters);
#endif // DUMP_DENSE_BLOCKS

    //    dumpAhmedMblockArray<ResultType>(blocks, blockCount);

    // // Dump timing data of individual chunks
    //    std::cout << "\nChunks:\n";
    //    for (int i = 0; i < leafClusterCount; ++i)
    //        if (chunkStats[i].valid) {
    //            int blockIndex = leafClusters[i]->getidx();
    //            std::cout << chunkStats[i].chunkStart << "\t"
    //                      << chunkStats[i].chunkSize << "\t"
    //                      << (chunkStats[i].startTime - loopStart).seconds() << "\t"
    //                      << (chunkStats[i].endTime - loopStart).seconds() << "\t"
    //                      << (chunkStats[i].endTime - chunkStats[i].startTime).seconds() << "\t"
    //                      << blocks[blockIndex]->getn1() << "\t"
    //                      << blocks[blockIndex]->getn2() << "\t"
    //                      << blocks[blockIndex]->islwr() << "\t"
    //                      << (blocks[blockIndex]->islwr() ? blocks[blockIndex]->rank() : 0) << "\n";
    //        }

    {
        size_t origMemory = sizeof(ResultType) * testDofCount * trialDofCount;
        size_t ahmedMemory = sizeH(bemBlclusterTree.get(), blocks.get());
        int maximumRank = Hmax_rank(bemBlclusterTree.get(), blocks.get());
        if (verbosityAtLeastDefault)
            std::cout << "\nNeeded storage: "
                      << ahmedMemory / 1024. / 1024. << " MB.\n"
                      << "Without approximation: "
                      << origMemory / 1024. / 1024. << " MB.\n"
                      << "Compressed to "
                      << (100. * ahmedMemory) / origMemory << "%.\n"
                      << "Maximum rank: " << maximumRank << ".\n"
                      << std::endl;

        if (acaOptions.outputPostscript) {
            if (verbosityAtLeastDefault)
                std::cout << "Writing matrix partition ..." << std::flush;
            std::ofstream os(acaOptions.outputFname.c_str());
            if (symmetric)
                // psoutputHeH() seems to work also for symmetric matrices
                psoutputHeH(os, bemBlclusterTree.get(),
                            trialDofCount, blocks.get());
            else
                psoutputGeH(os, bemBlclusterTree.get(),
                            std::max(testDofCount, trialDofCount), blocks.get());
            os.close();
            if (verbosityAtLeastDefault)
                std::cout << " done." << std::endl;
        }
    }

    int outSymmetry = NO_SYMMETRY;
    if (symmetric) {
        outSymmetry = SYMMETRIC;
        if (!boost::is_complex<ResultType>())
            outSymmetry |= HERMITIAN;
    }
    typedef DiscreteAcaBoundaryOperator<ResultType> DiscreteAcaLinOp;
    std::auto_ptr<DiscreteAcaLinOp> acaOp(
                new DiscreteAcaLinOp(testDofCount, trialDofCount,
                                     acaOptions.eps,
                                     acaOptions.maximumRank,
                                     outSymmetry,
                                     bemBlclusterTree, blocks,
                                     *trial_o2pPermutation, // domain
                                     *test_o2pPermutation, // range
                                     parallelOptions));
    return acaOp;
}

#endif

} // namespace

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
AcaGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        const std::vector<LocalAssemblerForBoundaryOperators*>& localAssemblers,
        const std::vector<const DiscreteBndOp*>& sparseTermsToAdd,
        const std::vector<ResultType>& denseTermMultipliers,
        const std::vector<ResultType>& sparseTermMultipliers,
        const AssemblyOptions& options,
        int symmetry)
{
#ifdef WITH_AHMED
    typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
    typedef ExtendedBemCluster<AhmedDofType> AhmedBemCluster;
    typedef bemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;
    typedef DiscreteAcaBoundaryOperator<ResultType> DiscreteAcaLinOp;

    const AcaOptions& acaOptions = options.acaOptions();
    const bool indexWithGlobalDofs = acaOptions.globalAssemblyBeforeCompression;
    const bool verbosityAtLeastDefault =
            (options.verbosityLevel() >= VerbosityLevel::DEFAULT);
    const bool verbosityAtLeastHigh =
            (options.verbosityLevel() >= VerbosityLevel::HIGH);

    // Currently we don't support Hermitian ACA operators. This is because we
    // don't have the means to really test them -- we would need complex-valued
    // basis functions for that. (Assembly of such a matrix would be very easy
    // -- just change complex_sym from true to false in the call to apprx_sym()
    // in AcaWeakFormAssemblerLoopBody::operator() -- but operations on
    // symmetric/Hermitian matrices are not always trivial and we do need to be
    // able to test them properly.)
    bool symmetric = symmetry & SYMMETRIC;
    if (symmetry & HERMITIAN && !(symmetry & SYMMETRIC) &&
            verbosityAtLeastDefault)
        std::cout << "Warning: assembly of non-symmetric Hermitian H-matrices "
                     "is not supported yet. A general H-matrix will be assembled"
                  << std::endl;

#ifndef WITH_TRILINOS
    if (!indexWithGlobalDofs)
        throw std::runtime_error("AcaGlobalAssembler::assembleDetachedWeakForm(): "
                                 "ACA assembly with globalAssemblyBeforeCompression "
                                 "set to false requires BEM++ to be linked with "
                                 "Trilinos");
#endif // WITH_TRILINOS

    const size_t testDofCount = indexWithGlobalDofs ?
                testSpace.globalDofCount() : testSpace.flatLocalDofCount();
    const size_t trialDofCount = indexWithGlobalDofs ?
                trialSpace.globalDofCount() : trialSpace.flatLocalDofCount();

    if (symmetric && testDofCount != trialDofCount)
        throw std::invalid_argument("AcaGlobalAssembler::assembleDetachedWeakForm(): "
                                    "you cannot generate a symmetric weak form "
                                    "using test and trial spaces with different "
                                    "numbers of DOFs");

    // o2p: map of original indices to permuted indices
    // p2o: map of permuted indices to original indices
    typedef ClusterConstructionHelper<BasisFunctionType> CCH;
    shared_ptr<AhmedBemCluster> testClusterTree;
    shared_ptr<IndexPermutation> test_o2pPermutation, test_p2oPermutation;
    CCH::constructBemCluster(testSpace, indexWithGlobalDofs, acaOptions,
                             testClusterTree,
                             test_o2pPermutation, test_p2oPermutation);
    shared_ptr<AhmedBemCluster> trialClusterTree;
    shared_ptr<IndexPermutation> trial_o2pPermutation, trial_p2oPermutation;
    if (symmetric || &testSpace == &trialSpace) {
        trialClusterTree = testClusterTree;
        trial_o2pPermutation = test_o2pPermutation;
        trial_p2oPermutation = test_p2oPermutation;
    } else
        CCH::constructBemCluster(trialSpace, indexWithGlobalDofs, acaOptions,
                                 trialClusterTree,
                                 trial_o2pPermutation, trial_p2oPermutation);

//    // Export VTK plots showing the disctribution of leaf cluster ids
//    std::vector<unsigned int> testClusterIds;
//    getClusterIds(*testClusterTree, test_p2oPermutation->permutedIndices(), testClusterIds);
//    testSpace.dumpClusterIds("testClusterIds", testClusterIds,
//                             indexWithGlobalDofs ? GLOBAL_DOFS : FLAT_LOCAL_DOFS);
//    std::vector<unsigned int> trialClusterIds;
//    getClusterIds(*trialClusterTree, trial_p2oPermutation->permutedIndices(), trialClusterIds);
//    trialSpace.dumpClusterIds("trialClusterIds", trialClusterIds,
//                              indexWithGlobalDofs ? GLOBAL_DOFS : FLAT_LOCAL_DOFS);

    if (verbosityAtLeastHigh)
        std::cout << "Test cluster count: " << testClusterTree->getncl()
                  << "\nTrial cluster count: " << trialClusterTree->getncl()
                  << std::endl;

    unsigned int blockCount = 0;
    shared_ptr<AhmedBemBlcluster> bemBlclusterTree(
                CCH::constructBemBlockCluster(acaOptions, symmetric,
                                              *testClusterTree, *trialClusterTree,
                                              blockCount).release());

    if (verbosityAtLeastHigh)
        std::cout << "Mblock count: " << blockCount << std::endl;

    std::vector<unsigned int> p2oTestDofs =
        test_p2oPermutation->permutedIndices();
    std::vector<unsigned int> p2oTrialDofs =
        trial_p2oPermutation->permutedIndices();
    assert(p2oTestDofs.size() == testDofCount);
    assert(p2oTrialDofs.size() == trialDofCount);

#ifdef DUMP_DENSE_BLOCKS
    std::vector<Point3D<CoordinateType> > testDofCenters, trialDofCenters;
    if (indexWithGlobalDofs) {
        testSpace.getGlobalDofPositions(testDofCenters);
        trialSpace.getGlobalDofPositions(trialDofCenters);
    } else {
        testSpace.getFlatLocalDofPositions(testDofCenters);
        trialSpace.getFlatLocalDofPositions(trialDofCenters);
    }
#endif // DUMP_DENSE_BLOCKS

    typedef WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>
            AcaAssemblyHelper;
    // TODO: It might be better (more efficient and elegant)
    // to pass p2oPermutation than p2oDofs.
    // Also, it might be more logical to rename IndexPermutation to IndexMapping
    // and permute/unpermute to map/unmap.
    AcaAssemblyHelper helper(
                testSpace, trialSpace, p2oTestDofs, p2oTrialDofs,
                localAssemblers, sparseTermsToAdd,
                denseTermMultipliers, sparseTermMultipliers, options);

    std::auto_ptr<DiscreteAcaBoundaryOperator<ResultType> > acaOp =
    assembleAcaOperator<AcaAssemblyHelper, BasisFunctionType, ResultType>(
                helper, bemBlclusterTree,
                options.parallelizationOptions(), options.acaOptions(),
                verbosityAtLeastDefault, symmetric,
                test_o2pPermutation, trial_o2pPermutation
#ifdef DUMP_DENSE_BLOCKS
                ,
                test_p2oPermutation, trial_p2oPermutation,
                testDofCenters, trialDofCenters
#endif // DUMP_DENSE_BLOCKS
                );

    std::auto_ptr<DiscreteBndOp> result;
    if (indexWithGlobalDofs)
        result = acaOp;
    else {
#ifdef WITH_TRILINOS
        // without Trilinos, this code will never be reached -- an exception
        // will be thrown earlier in this function
        typedef DiscreteBoundaryOperatorComposition<ResultType> DiscreteBndOpComp;
        shared_ptr<DiscreteBndOp> acaOpShared(acaOp.release());
        shared_ptr<DiscreteBndOp> trialGlobalToLocal =
                constructOperatorMappingGlobalToFlatLocalDofs<
                BasisFunctionType, ResultType>(trialSpace);
        shared_ptr<DiscreteBndOp> testLocalToGlobal =
                constructOperatorMappingFlatLocalToGlobalDofs<
                BasisFunctionType, ResultType>(testSpace);
        shared_ptr<DiscreteBndOp> tmp(
                    new DiscreteBndOpComp(acaOpShared, trialGlobalToLocal));
        result.reset(new DiscreteBndOpComp(testLocalToGlobal, tmp));
#endif // WITH_TRILINOS
    }
    return result;

#else // without Ahmed
    throw std::runtime_error("AcaGlobalAssembler::assembleDetachedWeakForm(): "
                             "To enable assembly in ACA mode, recompile BEM++ "
                             "with the symbol WITH_AHMED defined.");
#endif // WITH_AHMED
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
AcaGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        LocalAssemblerForBoundaryOperators& localAssembler,
        const AssemblyOptions& options,
        int symmetry)
{
    std::vector<LocalAssemblerForBoundaryOperators*> localAssemblers(
                1, &localAssembler);
    std::vector<const DiscreteBndOp*> sparseTermsToAdd;
    std::vector<ResultType> denseTermsMultipliers(1, 1.0);
    std::vector<ResultType> sparseTermsMultipliers;

    return assembleDetachedWeakForm(testSpace, trialSpace, localAssemblers,
                            sparseTermsToAdd,
                            denseTermsMultipliers,
                            sparseTermsMultipliers,
                            options, symmetry);
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
AcaGlobalAssembler<BasisFunctionType, ResultType>::assemblePotentialOperator(
        const arma::Mat<CoordinateType>& points,
        const Space<BasisFunctionType>& trialSpace,
        const std::vector<LocalAssemblerForPotentialOperators*>& localAssemblers,
        const std::vector<ResultType>& termMultipliers,
        const EvaluationOptions& options)
{
    const int symmetric = false;

#ifdef WITH_AHMED
    typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
    typedef ExtendedBemCluster<AhmedDofType> AhmedBemCluster;
    typedef bemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;
    typedef DiscreteAcaBoundaryOperator<ResultType> DiscreteAcaLinOp;

    const AcaOptions& acaOptions = options.acaOptions();
    const bool indexWithGlobalDofs = acaOptions.globalAssemblyBeforeCompression;
    const bool verbosityAtLeastDefault =
            (options.verbosityLevel() >= VerbosityLevel::DEFAULT);
    const bool verbosityAtLeastHigh =
            (options.verbosityLevel() >= VerbosityLevel::HIGH);

#ifndef WITH_TRILINOS
    if (!indexWithGlobalDofs)
        throw std::runtime_error("AcaGlobalAssembler::assemblePotentialOperator(): "
                                 "ACA assembly with globalAssemblyBeforeCompression "
                                 "set to false requires BEM++ to be linked with "
                                 "Trilinos");
#endif // WITH_TRILINOS

    if (localAssemblers.empty())
        throw std::runtime_error("AcaGlobalAssembler::assemblePotentialOperator(): "
                                 "the 'localAssemblers' vector must not be empty");

    const size_t pointCount = points.n_cols;
    const int componentCount = localAssemblers[0]->resultDimension();
    const size_t testDofCount = pointCount * componentCount;
    const size_t trialDofCount = indexWithGlobalDofs ?
                trialSpace.globalDofCount() : trialSpace.flatLocalDofCount();

    // o2p: map of original indices to permuted indices
    // p2o: map of permuted indices to original indices
    typedef ClusterConstructionHelper<BasisFunctionType> CCH;
    shared_ptr<AhmedBemCluster> testClusterTree;
    shared_ptr<IndexPermutation> test_o2pPermutation, test_p2oPermutation;
    CCH::constructBemCluster(points, componentCount, acaOptions,
                             testClusterTree,
                             test_o2pPermutation, test_p2oPermutation);
    shared_ptr<AhmedBemCluster> trialClusterTree;
    shared_ptr<IndexPermutation> trial_o2pPermutation, trial_p2oPermutation;
    CCH::constructBemCluster(trialSpace, indexWithGlobalDofs, acaOptions,
                             trialClusterTree,
                             trial_o2pPermutation, trial_p2oPermutation);

    // Print the distribution of cluster ids
#ifdef DUMP_DENSE_BLOCKS
    std::vector<Point3D<CoordinateType> > testDofCenters, trialDofCenters;
    CCH::getComponentDofPositions(points, componentCount, testDofCenters);
    if (indexWithGlobalDofs)
        trialSpace.getGlobalDofPositions(trialDofCenters);
    else
        trialSpace.getFlatLocalDofPositions(trialDofCenters);
#endif // DUMP_DENSE_BLOCKS

    if (verbosityAtLeastHigh)
        std::cout << "Test cluster count: " << testClusterTree->getncl()
                  << "\nTrial cluster count: " << trialClusterTree->getncl()
                  << std::endl;

    unsigned int blockCount = 0;
    shared_ptr<AhmedBemBlcluster> bemBlclusterTree(
                CCH::constructBemBlockCluster(acaOptions, false /* symmetric */,
                                              *testClusterTree, *trialClusterTree,
                                              blockCount).release());

    if (verbosityAtLeastHigh)
        std::cout << "Mblock count: " << blockCount << std::endl;

    std::vector<unsigned int> p2oPoints =
        test_p2oPermutation->permutedIndices();
    std::vector<unsigned int> p2oTrialDofs =
        trial_p2oPermutation->permutedIndices();
    typedef PotentialOperatorAcaAssemblyHelper<BasisFunctionType, ResultType>
            AcaAssemblyHelper;
    AcaAssemblyHelper helper(points, trialSpace, p2oPoints, p2oTrialDofs,
                             localAssemblers, termMultipliers, options);

    std::auto_ptr<DiscreteAcaBoundaryOperator<ResultType> > acaOp =
    assembleAcaOperator<AcaAssemblyHelper, BasisFunctionType, ResultType>(
                helper, bemBlclusterTree,
                options.parallelizationOptions(), options.acaOptions(),
                verbosityAtLeastDefault, symmetric,
                test_o2pPermutation, trial_o2pPermutation
#ifdef DUMP_DENSE_BLOCKS
                ,
                test_p2oPermutation, trial_p2oPermutation,
                testDofCenters, trialDofCenters
#endif // DUMP_DENSE_BLOCKS
                );

    std::auto_ptr<DiscreteBndOp> result;
    if (indexWithGlobalDofs)
        result = acaOp;
    else {
#ifdef WITH_TRILINOS
        // without Trilinos, this code will never be reached -- an exception
        // will be thrown earlier in this function
        typedef DiscreteBoundaryOperatorComposition<ResultType> DiscreteBndOpComp;
        shared_ptr<DiscreteBndOp> acaOpShared(acaOp.release());
        shared_ptr<DiscreteBndOp> trialGlobalToLocal =
                constructOperatorMappingGlobalToFlatLocalDofs<
                BasisFunctionType, ResultType>(trialSpace);
        result.reset(new DiscreteBndOpComp(acaOpShared, trialGlobalToLocal));
#endif // WITH_TRILINOS
    }
    return result;

#else // without Ahmed
    throw std::runtime_error("AcaGlobalAssembler::assemblePotentialOperator(): "
                             "To enable assembly in ACA mode, recompile BEM++ "
                             "with the symbol WITH_AHMED defined.");
#endif // WITH_AHMED
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
AcaGlobalAssembler<BasisFunctionType, ResultType>::assemblePotentialOperator(
        const arma::Mat<CoordinateType>& points,
        const Space<BasisFunctionType>& trialSpace,
        LocalAssemblerForPotentialOperators& localAssembler,
        const EvaluationOptions& options)
{
    std::vector<LocalAssemblerForPotentialOperators*> localAssemblers(
                1, &localAssembler);
    std::vector<ResultType> termMultipliers(1, 1.0);

    return assemblePotentialOperator(points, trialSpace, localAssemblers,
                                     termMultipliers, options);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AcaGlobalAssembler);

} // namespace Bempp
