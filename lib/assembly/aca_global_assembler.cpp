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
#include "index_permutation.hpp"
#include "discrete_boundary_operator_composition.hpp"
#include "discrete_sparse_boundary_operator.hpp"

#include "../common/armadillo_fwd.hpp"
#include "../common/auto_timer.hpp"
#include "../common/boost_shared_array_fwd.hpp"
#include "../common/chunk_statistics.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_operators.hpp"
#include "../fiber/serial_blas_region.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../space/space.hpp"

#include <stdexcept>
#include <iostream>

#include <boost/type_traits/is_complex.hpp>

#include <tbb/atomic.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_queue.h>

#ifdef WITH_AHMED
#include "ahmed_aux.hpp"
#include "discrete_aca_boundary_operator.hpp"
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
            bool verbose,
            LeafClusterIndexQueue& leafClusterIndexQueue,
            bool hermitian,
            std::vector<ChunkStatistics>& stats) :
        m_helper(helper),
        m_leafClusters(leafClusters), m_blocks(blocks),
        m_options(options), m_done(done), m_verbose(verbose),
        m_leafClusterIndexQueue(leafClusterIndexQueue),
        m_hermitian(hermitian),
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
            if (m_hermitian)
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
            if (m_verbose)
                progressbar(std::cout, TEXT, (++m_done) - 1,
                            m_leafClusters.size(), HASH_COUNT, true);
        }

    }

private:
    WeakFormAcaAssemblyHelper<BasisFunctionType, ResultType>& m_helper;
    AhmedLeafClusterArray& m_leafClusters;
    boost::shared_array<AhmedMblock*> m_blocks;
    const AcaOptions& m_options;
    tbb::atomic<size_t>& m_done;
    bool m_verbose;
    LeafClusterIndexQueue& m_leafClusterIndexQueue;
    bool m_hermitian;
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
#endif

} // namespace

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
AcaGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
        const Space<BasisFunctionType>& testSpace,
        const Space<BasisFunctionType>& trialSpace,
        const std::vector<LocalAssembler*>& localAssemblers,
        const std::vector<const DiscreteBndOp*>& sparseTermsToAdd,
        const std::vector<ResultType>& denseTermsMultipliers,
        const std::vector<ResultType>& sparseTermsMultipliers,
        const AssemblyOptions& options,
        bool hermitian)
{
#ifdef WITH_AHMED
    typedef AhmedDofWrapper<CoordinateType> AhmedDofType;
    typedef DiscreteAcaBoundaryOperator<ResultType> DiscreteAcaLinOp;

    const AcaOptions& acaOptions = options.acaOptions();
    const bool indexWithGlobalDofs = acaOptions.globalAssemblyBeforeCompression;
    const bool verbosityDefault =
            (options.verbosityLevel() >= VerbosityLevel::DEFAULT);
    const bool verbosityHigh =
            (options.verbosityLevel() >= VerbosityLevel::HIGH);

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

    if (hermitian && testDofCount != trialDofCount)
        throw std::invalid_argument("AcaGlobalAssembler::assembleDetachedWeakForm(): "
                                    "you cannot generate a Hermitian weak form "
                                    "using test and trial spaces with different "
                                    "numbers of DOFs");

    // o2p: map of original indices to permuted indices
    // p2o: map of permuted indices to original indices
    std::vector<unsigned int> o2pTestDofs(testDofCount);
    std::vector<unsigned int> p2oTestDofs(testDofCount);
    std::vector<unsigned int> o2pTrialDofs(trialDofCount);
    std::vector<unsigned int> p2oTrialDofs(trialDofCount);
    for (size_t i = 0; i < testDofCount; ++i)
        o2pTestDofs[i] = i;
    for (size_t i = 0; i < testDofCount; ++i)
        p2oTestDofs[i] = i;
    for (size_t i = 0; i < trialDofCount; ++i)
        o2pTrialDofs[i] = i;
    for (size_t i = 0; i < trialDofCount; ++i)
        p2oTrialDofs[i] = i;

    std::vector<Point3D<CoordinateType> > trialDofCenters, testDofCenters;
    if (indexWithGlobalDofs) {
        trialSpace.getGlobalDofPositions(trialDofCenters);
        testSpace.getGlobalDofPositions(testDofCenters);
    } else {
        trialSpace.getFlatLocalDofPositions(trialDofCenters);
        testSpace.getFlatLocalDofPositions(testDofCenters);
    }

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

    if (verbosityHigh)
        std::cout << "Test cluster count: " << testClusterTree.getncl()
                  << "\nTrial cluster count: " << trialClusterTree.getncl()
                  << std::endl;
    //    std::cout << "o2pTest:\n" << o2pTestDofs << std::endl;
    //    std::cout << "p2oTest:\n" << p2oTestDofs << std::endl;

    typedef bemblcluster<AhmedDofType, AhmedDofType> AhmedBemBlcluster;
    std::auto_ptr<AhmedBemBlcluster> bemBlclusterTree(
                new AhmedBemBlcluster(0, 0, testDofCount, trialDofCount));
    unsigned int blockCount = 0;
    if (hermitian)
        bemBlclusterTree->subdivide_sym(&testClusterTree,
                                        acaOptions.eta * acaOptions.eta,
                                        blockCount);
    else
        bemBlclusterTree->subdivide(&testClusterTree, &trialClusterTree,
                                    acaOptions.eta * acaOptions.eta,
                                    blockCount);

    if (verbosityHigh)
        std::cout << "Double cluster count: " << blockCount << std::endl;

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

    const ParallelizationOptions& parallelOptions =
            options.parallelizationOptions();
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

    if (verbosityDefault)
        std::cout << "About to start the ACA assembly loop" << std::endl;
    tbb::tick_count loopStart = tbb::tick_count::now();    
    {
        Fiber::SerialBlasRegion region; // if possible, ensure that BLAS is single-threaded
        tbb::parallel_for(tbb::blocked_range<size_t>(0, leafClusterCount),
                          Body(helper, leafClusters, blocks, acaOptions, done,
                               verbosityDefault,
                               leafClusterIndexQueue, hermitian, chunkStats));
    }
    tbb::tick_count loopEnd = tbb::tick_count::now();
    if (verbosityDefault) {
        std::cout << "\n"; // the progress bar doesn't print the final \n
        std::cout << "ACA loop took " << (loopEnd - loopStart).seconds() << " s"
                  << std::endl;
    }

    // TODO: parallelise!
    if (acaOptions.recompress) {
        if (verbosityDefault)
            std::cout << "About to start ACA agglomeration" << std::endl;
        agglH(bemBlclusterTree.get(), blocks.get(),
              acaOptions.eps, acaOptions.maximumRank);
        if (verbosityDefault)
            std::cout << "Agglomeration finished" << std::endl;
    }

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
        if (verbosityDefault)
            std::cout << "\nNeeded storage: "
                      << ahmedMemory / 1024. / 1024. << " MB.\n"
                      << "Without approximation: "
                      << origMemory / 1024. / 1024. << " MB.\n"
                      << "Compressed to "
                      << (100. * ahmedMemory) / origMemory << "%.\n"
                      << std::endl;

        if (acaOptions.outputPostscript) {
            if (verbosityDefault)
                std::cout << "Writing matrix partition ..." << std::flush;
            std::ofstream os(acaOptions.outputFname.c_str());
            if (hermitian)
#ifdef AHMED_PRERELEASE
                psoutputHSym(os, bemBlclusterTree.get(), testDofCount, blocks.get());
#else
                psoutputHeH(os, bemBlclusterTree.get(), testDofCount, blocks.get());
#endif
            else
#ifdef AHMED_PRERELEASE
                psoutputH(os, bemBlclusterTree.get(), testDofCount, blocks.get());
#else
                psoutputGeH(os, bemBlclusterTree.get(), testDofCount, blocks.get());
#endif
            os.close();
            if (verbosityDefault)
                std::cout << " done." << std::endl;
        }
    }

    unsigned int symmetry = NO_SYMMETRY;
    if (hermitian) {
        symmetry |= HERMITIAN;
        if (boost::is_complex<ResultType>())
            symmetry |= SYMMETRIC;
    }
    std::auto_ptr<DiscreteBndOp> acaOp(
                new DiscreteAcaLinOp(testDofCount, trialDofCount,
                                     acaOptions.maximumRank,
                                     Symmetry(symmetry),
                                     bemBlclusterTree, blocks,
                                     IndexPermutation(o2pTrialDofs),
                                     IndexPermutation(o2pTestDofs)));

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
        LocalAssembler& localAssembler,
        const AssemblyOptions& options,
        bool hermitian)
{
    std::vector<LocalAssembler*> localAssemblers(1, &localAssembler);
    std::vector<const DiscreteBndOp*> sparseTermsToAdd;
    std::vector<ResultType> denseTermsMultipliers(1, 1.0);
    std::vector<ResultType> sparseTermsMultipliers;

    return assembleDetachedWeakForm(testSpace, trialSpace, localAssemblers,
                            sparseTermsToAdd,
                            denseTermsMultipliers,
                            sparseTermsMultipliers,
                            options, hermitian);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AcaGlobalAssembler);

} // namespace Bempp
