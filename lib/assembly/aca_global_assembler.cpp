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

#include "config_ahmed.hpp"

#include "aca_global_assembler.hpp"

#include "assembly_options.hpp"
#include "index_permutation.hpp"

#include "../common/auto_timer.hpp"
#include "../common/chunk_statistics.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_operators.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../space/space.hpp"

#include <armadillo>
#include <boost/shared_array.hpp>
#include <stdexcept>
#include <iostream>

#include <tbb/atomic.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_queue.h>

#ifdef WITH_AHMED
#include "ahmed_aux.hpp"
#include "discrete_aca_linear_operator.hpp"
#include "scattered_range.hpp"
#include "weak_form_aca_assembly_helper.hpp"
#endif

namespace Bempp
{

// Body of parallel loop
namespace
{

#ifdef WITH_AHMED
template <typename BasisFunctionType, typename ResultType>
class AcaWeakFormAssemblerLoopBody
{
    typedef typename Fiber::ScalarTraits<ResultType>::RealType CoordinateType;
    typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
    typedef bemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;
    typedef mblock<typename AhmedTypeTraits<ResultType>::Type> AhmedMblock;
public:
    typedef tbb::concurrent_queue<size_t> LeafClusterIndexQueue;

    AcaWeakFormAssemblerLoopBody(
            WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>& helper,
            AhmedLeafClusterArray& leafClusters,
            boost::shared_array<AhmedMblock*> blocks,
            const AcaOptions& options,
            tbb::atomic<size_t>& done,
            LeafClusterIndexQueue& leafClusterIndexQueue,
            bool symmetric,
            std::vector<ChunkStatistics>& stats) :
        m_helper(helper),
        m_leafClusters(leafClusters), m_blocks(blocks),
        m_options(options), m_done(done),
        m_leafClusterIndexQueue(leafClusterIndexQueue),
        m_symmetric(symmetric),
        m_stats(stats)
    {
    }

    template <typename Range>
    void operator() (const Range& r) const {
        const char* TEXT = "Approximating ... ";
        for (typename Range::const_iterator i = r.begin(); i != r.end(); ++i) {
            size_t leafClusterIndex = -1;
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
#ifdef AHMED_PRERELEASE
                apprx_sym(m_helper, m_blocks[cluster->getidx()],
                          cluster, m_options.eps, m_options.maximumRank);
#else
                apprx_sym(m_helper, m_blocks[cluster->getidx()],
                          cluster, m_options.eps, m_options.maximumRank,
                          true /* complex_sym */);
#endif
            else
                apprx_unsym(m_helper, m_blocks[cluster->getidx()],
                            cluster, m_options.eps, m_options.maximumRank);
            m_stats[leafClusterIndex].endTime = tbb::tick_count::now();
            // TODO: recompress
            const int HASH_COUNT = 20;
            progressbar(std::cout, TEXT, (++m_done) - 1,
                        m_leafClusters.size(), HASH_COUNT, true);
        }

    }

private:
    mutable WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>& m_helper;
    AhmedLeafClusterArray& m_leafClusters;
    boost::shared_array<AhmedMblock*> m_blocks;
    const AcaOptions& m_options;
    mutable tbb::atomic<size_t>& m_done;
    mutable LeafClusterIndexQueue& m_leafClusterIndexQueue;
    bool m_symmetric;
    mutable std::vector<ChunkStatistics>& m_stats;
};

void reallyGetClusterIds(const cluster& clusterTree,
                         const std::vector<unsigned int>& p2oDofs,
                         std::vector<unsigned int>& clusterIds,
                         unsigned int& id)
{
    if (clusterTree.isleaf())
        for (int nDof = clusterTree.getnbeg(); nDof < clusterTree.getnend(); ++nDof)
            clusterIds[p2oDofs[nDof]] = id;
    else
        for (int nSon = 0; nSon < clusterTree.getns(); ++nSon)
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
#endif

} // namespace

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
AcaGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        const std::vector<LocalAssembler*>& localAssemblers,
        const std::vector<const DiscreteLinOp*>& sparseTermsToAdd,
        const std::vector<ResultType>& denseTermsMultipliers,
        const std::vector<ResultType>& sparseTermsMultipliers,
        const AssemblyOptions& options,
        bool symmetric)
{
#ifdef WITH_AHMED
    typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
    typedef DiscreteAcaLinearOperator<ResultType> DiscreteAcaLinOp;

    const AcaOptions& acaOptions = options.acaOptions();

    const int testDofCount = testSpace.globalDofCount();
    const int trialDofCount = trialSpace.globalDofCount();

    if (symmetric && testDofCount != trialDofCount)
        throw std::invalid_argument("AcaGlobalAssembler::assembleDetachedWeakForm(): "
                                    "you cannot generate a symmetric weak form "
                                    "using test and trial spaces with different "
                                    "numbers of DOFs");

#ifndef NDEBUG
    std::cout << "Generating cluster trees... " << std::endl;
#endif
    // o2p: map of original indices to permuted indices
    // p2o: map of permuted indices to original indices
    std::vector<unsigned int> o2pTestDofs(testDofCount);
    std::vector<unsigned int> p2oTestDofs(testDofCount);
    std::vector<unsigned int> o2pTrialDofs(trialDofCount);
    std::vector<unsigned int> p2oTrialDofs(trialDofCount);
    for (unsigned int i = 0; i < testDofCount; ++i)
        o2pTestDofs[i] = i;
    for (unsigned int i = 0; i < testDofCount; ++i)
        p2oTestDofs[i] = i;
    for (unsigned int i = 0; i < trialDofCount; ++i)
        o2pTrialDofs[i] = i;
    for (unsigned int i = 0; i < trialDofCount; ++i)
        p2oTrialDofs[i] = i;

    std::vector<Point3D<CoordinateType> > trialDofCenters, testDofCenters;
    trialSpace.globalDofPositions(trialDofCenters);
    testSpace.globalDofPositions(testDofCenters);

    // Use static_cast to convert from a pointer to Point3D to a pointer to its
    // descendant AhmedDofWrapper, which does not contain any new data members,
    // but just one additional method (the two structs should therefore be
    // binary compatible)
    const AhmedDofType* ahmedTrialDofCenters =
            static_cast<AhmedDofType*>(&trialDofCenters[0]);
    const AhmedDofType* ahmedTestDofCenters =
            static_cast<AhmedDofType*>(&testDofCenters[0]);

    // NOTE: Ahmed uses names "op_perm" and "po_perm", which
    // correspond to BEM++'s "p2o" and "o2p", NOT the other way round.
    ExtendedBemCluster<const AhmedDofType> testClusterTree(
                ahmedTestDofCenters, &p2oTestDofs[0],
                0, testDofCount, acaOptions.maximumBlockSize);
    testClusterTree.createClusterTree(
                acaOptions.minimumBlockSize,
                &p2oTestDofs[0], &o2pTestDofs[0]);
    ExtendedBemCluster<const AhmedDofType> trialClusterTree(
                ahmedTrialDofCenters, &p2oTrialDofs[0],
                0, trialDofCount, acaOptions.maximumBlockSize);
    trialClusterTree.createClusterTree(
                acaOptions.minimumBlockSize,
                &p2oTrialDofs[0], &o2pTrialDofs[0]);

    // // Export VTK plots showing the disctribution of leaf cluster ids
    // std::vector<unsigned int> testClusterIds;
    // getClusterIds(testClusterTree, p2oTestDofs, testClusterIds);
    // testSpace.dumpClusterIds("testClusterIds", testClusterIds);
    // std::vector<unsigned int> trialClusterIds;
    // getClusterIds(trialClusterTree, p2oTrialDofs, trialClusterIds);
    // trialSpace.dumpClusterIds("trialClusterIds", trialClusterIds);

#ifndef NDEBUG
    std::cout << "Test cluster count: " << testClusterTree.getncl()
              << "\nTrial cluster count: " << trialClusterTree.getncl()
              << std::endl;
    //    std::cout << "o2pTest:\n" << o2pTestDofs << std::endl;
    //    std::cout << "p2oTest:\n" << p2oTestDofs << std::endl;
#endif

    typedef bemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;
    std::auto_ptr<AhmedBemBlcluster> bemBlclusterTree(
                new AhmedBemBlcluster(0, 0, testDofCount, trialDofCount));
    unsigned int blockCount = 0;
    if (symmetric)
        bemBlclusterTree->subdivide_sym(&testClusterTree,
                                        acaOptions.eta * acaOptions.eta,
                                        blockCount);
    else
        bemBlclusterTree->subdivide(&testClusterTree, &trialClusterTree,
                                    acaOptions.eta * acaOptions.eta,
                                    blockCount);

#ifndef NDEBUG
    std::cout << "Double cluster count: " << blockCount << std::endl;
#endif

    std::auto_ptr<DiscreteLinOp> result;

    WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>
            helper(testSpace, trialSpace, p2oTestDofs, p2oTrialDofs,
                   localAssemblers, sparseTermsToAdd,
                   denseTermsMultipliers, sparseTermsMultipliers, options);

    typedef mblock<typename AhmedTypeTraits<ResultType>::Type> AhmedMblock;
    boost::shared_array<AhmedMblock*> blocks =
            allocateAhmedMblockArray<ResultType>(bemBlclusterTree.get());

    // matgen_sqntl(helper, AhmedBemBlclusterTree.get(), AhmedBemBlclusterTree.get(),
    //              acaOptions.recompress, acaOptions.eps,
    //              acaOptions.maximumRank, blocks.get());

    // matgen_omp(helper, blockCount, AhmedBemBlclusterTree.get(),
    //            acaOptions.eps, acaOptions.maximumRank, blocks.get());

    // // Dump mblocks
    // const int mblockCount = AhmedBemBlclusterTree->nleaves();
    // for (int i = 0; i < mblockCount; ++i)
    //     if (blocks[i]->isdns())
    //     {
    //         char  buffer[1024];
    //         sprintf(buffer, "mblock-dns-%d-%d.txt",
    //                 blocks[i]->getn1(), blocks[i]->getn2());
    //         arma::Col<ResultType> block((ResultType*)blocks[i]->getdata(),
    //                                     blocks[i]->nvals());
    //         arma::diskio::save_raw_ascii(block, buffer);
    //     }
    //     else
    //     {
    //         char buffer[1024];
    //         sprintf(buffer, "mblock-lwr-%d-%d.txt",
    //                 blocks[i]->getn1(), blocks[i]->getn2());
    //         arma::Col<ResultType> block((ResultType*)blocks[i]->getdata(),
    //                                     blocks[i]->nvals());
    //         arma::diskio::save_raw_ascii(block, buffer);
    //     }

    AhmedLeafClusterArray leafClusters(bemBlclusterTree.get());
    leafClusters.sortAccordingToClusterSize();
    const size_t leafClusterCount = leafClusters.size();

    const ParallelisationOptions& parallelOptions =
            options.parallelisationOptions();
    int maxThreadCount = 1;
    if (parallelOptions.mode() == ParallelisationOptions::TBB)
    {
        if (parallelOptions.maxThreadCount() == ParallelisationOptions::AUTO)
            maxThreadCount = tbb::task_scheduler_init::automatic;
        else
            maxThreadCount = parallelOptions.maxThreadCount();
    }
    tbb::task_scheduler_init scheduler(maxThreadCount);
    tbb::atomic<size_t> done;
    done = 0;

    std::vector<ChunkStatistics> chunkStats(leafClusterCount);

    //    typedef AcaWeakFormAssemblerLoopBody<BasisFunctionType, ResultType> Body;
    //    // std::cout << "Loop start" << std::endl;
    //    tbb::tick_count loopStart = tbb::tick_count::now();
    // //    tbb::parallel_for(tbb::blocked_range<size_t>(0, leafClusterCount),
    // //                      Body(helper, leafClusters, blocks, acaOptions, done
    // //                           , chunkStats));
    //    tbb::parallel_for(ScatteredRange(0, leafClusterCount),
    //                      Body(helper, leafClusters, blocks, acaOptions, done
    //                           , chunkStats));
    //    tbb::tick_count loopEnd = tbb::tick_count::now();
    //    // std::cout << "Loop end" << std::endl;

    typedef AcaWeakFormAssemblerLoopBody<BasisFunctionType, ResultType> Body;
    typename Body::LeafClusterIndexQueue leafClusterIndexQueue;
    for (size_t i = 0; i < leafClusterCount; ++i)
        leafClusterIndexQueue.push(i);

    std::cout << "About to start the ACA assembly loop" << std::endl;
    tbb::tick_count loopStart = tbb::tick_count::now();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, leafClusterCount),
                      Body(helper, leafClusters, blocks, acaOptions, done,
                           leafClusterIndexQueue, symmetric, chunkStats));
    tbb::tick_count loopEnd = tbb::tick_count::now();

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

    std::cout << "\nTotal time: " << (loopEnd - loopStart).seconds() << std::endl;

    {
        size_t origMemory = sizeof(ResultType) * testDofCount * trialDofCount;
        size_t ahmedMemory = sizeH(bemBlclusterTree.get(), blocks.get());
        std::cout << "\nNeeded storage: " << ahmedMemory / 1024. / 1024. << " MB.\n"
                  << "Without approximation: " << origMemory / 1024. / 1024. << " MB.\n"
                  << "Compressed to " << (100. * ahmedMemory) / origMemory << "%.\n"
                  << std::endl;

        if (acaOptions.outputPostscript) {
            std::cout << "Writing matrix partition ..." << std::flush;
            std::ofstream os(acaOptions.outputFname.c_str());
            if (symmetric)
#if AHMED_PRERELEASE
                psoutputHSym(os, bemBlclusterTree.get(), testDofCount, blocks.get());
#else
                psoutputHeH(os, bemBlclusterTree.get(), testDofCount, blocks.get());
#endif
            else
#if AHMED_PRERELEASE
                psoutputH(os, bemBlclusterTree.get(), testDofCount, blocks.get());
#else
                psoutputGeH(os, bemBlclusterTree.get(), testDofCount, blocks.get());
#endif
            os.close();
            std::cout << " done." << std::endl;
        }
    }

    result = std::auto_ptr<DiscreteLinOp>(
                new DiscreteAcaLinOp(testDofCount, trialDofCount,
                                     acaOptions.maximumRank,
                                     symmetric,
                                     bemBlclusterTree, blocks,
                                     IndexPermutation(o2pTrialDofs),
                                     IndexPermutation(o2pTestDofs)));
    return result;
#else // without Ahmed
    throw std::runtime_error("AcaGlobalAssembler::assembleDetachedWeakForm(): "
                             "To enable assembly in ACA mode, recompile BEM++ "
                             "with the symbol WITH_AHMED defined.");
#endif
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
AcaGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        LocalAssembler& localAssembler,
        const AssemblyOptions& options,
        bool symmetric)
{
    std::vector<LocalAssembler*> localAssemblers(1, &localAssembler);
    std::vector<const DiscreteLinOp*> sparseTermsToAdd;
    std::vector<ResultType> denseTermsMultipliers(1, 1.0);
    std::vector<ResultType> sparseTermsMultipliers;

    return assembleDetachedWeakForm(testSpace, trialSpace, localAssemblers,
                            sparseTermsToAdd,
                            denseTermsMultipliers,
                            sparseTermsMultipliers,
                            options, symmetric);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AcaGlobalAssembler);

} // namespace Bempp
