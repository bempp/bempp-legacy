#include "../common/config_ahmed.hpp"

#include "aca_global_assembler.hpp"

#include "assembly_options.hpp"
#include "index_permutation.hpp"

#include "../common/auto_timer.hpp"
#include "../fiber/local_assembler_for_operators.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../space/space.hpp"

#include <armadillo>
#include <boost/shared_array.hpp>
#include <stdexcept>
#include <iostream>

#include <tbb/atomic.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

#ifdef WITH_AHMED
#include "ahmed_aux.hpp"
#include "discrete_aca_scalar_valued_linear_operator.hpp"
#include "weak_form_aca_assembly_helper.hpp"
#endif

namespace Bempp
{

// Body of parallel loop
namespace
{

#ifdef WITH_AHMED
template <typename ValueType>
class AcaWeakFormAssemblerLoopBody
{
public:
    typedef AhmedDofWrapper<ValueType> AhmedDofType;
    typedef bemblcluster<AhmedDofType, AhmedDofType> DoubleCluster;

    AcaWeakFormAssemblerLoopBody(
            WeakFormAcaAssemblyHelper<ValueType>& helper,
            AhmedLeafClusterArray& leafClusters,
            boost::shared_array<mblock<ValueType>*> blocks,
            const AcaOptions& options,
            tbb::atomic<size_t>& done) :
        m_helper(helper),
        m_leafClusters(leafClusters), m_blocks(blocks),
        m_options(options), m_done(done)
    {
    }

    void operator() (const tbb::blocked_range<size_t>& r) const {
        const char* TEXT = "Approximating ... ";
        for (size_t i = r.begin(); i != r.end(); ++i)
        {
            DoubleCluster* cluster = dynamic_cast<DoubleCluster*>(m_leafClusters[i]);
            apprx_unsym(m_helper, m_blocks[cluster->getidx()],
                        cluster, m_options.eps, m_options.maximumRank);
            // TODO: recompress
            const int HASH_COUNT = 20;
            progressbar(std::cout, TEXT, (++m_done) - 1,
                        m_leafClusters.size(), HASH_COUNT, true);
        }
    }

private:
    mutable WeakFormAcaAssemblyHelper<ValueType>& m_helper;
    size_t m_leafClusterCount;
    AhmedLeafClusterArray& m_leafClusters;
    boost::shared_array<mblock<ValueType>*> m_blocks;
    const AcaOptions& m_options;
    mutable tbb::atomic<size_t>& m_done;
};
#endif

} // namespace

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
AcaGlobalAssembler<ValueType>::assembleWeakForm(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        const std::vector<LocalAssembler*>& localAssemblers,
        const std::vector<const DiscreteLinOp*>& sparseTermsToAdd,
        const AssemblyOptions& options)
{
#ifdef WITH_AHMED
    typedef AhmedDofWrapper<ValueType> AhmedDofType;
    typedef DiscreteAcaScalarValuedLinearOperator<ValueType> DiscreteAcaLinOp;

    const AcaOptions& acaOptions = options.acaOptions();

    // Get the grid's leaf view so that we can iterate over elements
    std::auto_ptr<GridView> view = trialSpace.grid().leafView();

    // const int elementCount = view.entityCount(0);
    const int trialDofCount = trialSpace.globalDofCount();
    const int testDofCount = testSpace.globalDofCount();

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

    std::vector<Point3D<ValueType> > trialDofCenters, testDofCenters;
    trialSpace.globalDofPositions(trialDofCenters);
    testSpace.globalDofPositions(testDofCenters);

    // Use static_cast to convert from a pointer to Point3D to a pointer to its
    // descendant AhmedDofWrapper, which does not contain any new data members,
    // but just one additional method (the two structs should therefore be
    // binary compatible)
    const AhmedDofType* ahmedTrialDofCenters =
            static_cast<AhmedDofType*>(&trialDofCenters[0]);
    const AhmedDofType* ahmedTestDofCenters =
            static_cast<AhmedDofType*>(&trialDofCenters[0]);

    bemcluster<const AhmedDofType> testClusterTree(
                ahmedTestDofCenters, &o2pTestDofs[0],
                0, testDofCount);
    testClusterTree.createClusterTree(
                acaOptions.minimumBlockSize,
                &o2pTestDofs[0], &p2oTestDofs[0]);
    bemcluster<const AhmedDofType> trialClusterTree(
                ahmedTrialDofCenters, &o2pTrialDofs[0],
                0, trialDofCount);
    trialClusterTree.createClusterTree(
                acaOptions.minimumBlockSize,
                &o2pTrialDofs[0], &p2oTrialDofs[0]);

#ifndef NDEBUG
    std::cout << "Test cluster count: " << testClusterTree.getncl()
              << "\nTrial cluster count: " << trialClusterTree.getncl()
              << std::endl;
//    std::cout << "o2pTest:\n" << o2pTestDofs << std::endl;
//    std::cout << "p2oTest:\n" << p2oTestDofs << std::endl;
#endif

    typedef bemblcluster<AhmedDofType, AhmedDofType> DoubleCluster;
    std::auto_ptr<DoubleCluster> doubleClusterTree(
                new DoubleCluster(0, 0, testDofCount, trialDofCount));
    unsigned int blockCount = 0;
    doubleClusterTree->subdivide(&testClusterTree, &trialClusterTree,
                                 acaOptions.eta * acaOptions.eta,
                                 blockCount);

#ifndef NDEBUG
    std::cout << "Double cluster count: " << blockCount << std::endl;
#endif

    std::auto_ptr<DiscreteLinOp> result;

    WeakFormAcaAssemblyHelper<ValueType>
            helper(testSpace, trialSpace, p2oTestDofs, p2oTrialDofs,
                   localAssemblers, sparseTermsToAdd, options);

    boost::shared_array<mblock<ValueType>*> blocks =
            allocateAhmedMblockArray<ValueType>(doubleClusterTree.get());

//    matgen_sqntl(helper, doubleClusterTree.get(), doubleClusterTree.get(),
//                 acaOptions.recompress, acaOptions.eps,
//                 acaOptions.maximumRank, blocks);

//    matgen_omp(helper, blockCount, doubleClusterTree.get(),
//                   acaOptions.eps, acaOptions.maximumRank, blocks);

    AhmedLeafClusterArray leafClusters(doubleClusterTree.get());
    const size_t leafClusterCount = leafClusters.size();

    int maxThreadCount = 1;
    if (options.parallelism() == AssemblyOptions::TBB)
    {
        if (options.maxThreadCount() == AssemblyOptions::AUTO)
            maxThreadCount = tbb::task_scheduler_init::automatic;
        else
            maxThreadCount = options.maxThreadCount();
    }
    tbb::task_scheduler_init scheduler(maxThreadCount);
    tbb::atomic<size_t> done;
    done = 0;

    typedef AcaWeakFormAssemblerLoopBody<ValueType> Body;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, leafClusterCount),
                      Body(helper, leafClusters, blocks, acaOptions, done));

    {
        size_t origMemory = sizeof(ValueType) * testDofCount * trialDofCount;
        size_t ahmedMemory = sizeH(doubleClusterTree.get(), blocks.get());
        std::cout << "\nNeeded storage: " << ahmedMemory / 1024. / 1024. << " MB.\n"
                  << "Without approximation: " << origMemory / 1024. / 1024. << " MB.\n"
                  << "Compressed to " << (100. * ahmedMemory) / origMemory << "%.\n"
                  << std::endl;

        std::cout << "Writing matrix partition ..." << std::flush;
        std::ofstream os("aca.ps");
        psoutputH(os, doubleClusterTree.get(), testDofCount, blocks.get());
        os.close();
        std::cout << " done." << std::endl;
    }

    result = std::auto_ptr<DiscreteLinOp>(
                new DiscreteAcaLinOp(testDofCount, trialDofCount,
                                     acaOptions.maximumRank,
                                     doubleClusterTree, blocks,
                                     IndexPermutation(o2pTestDofs),
                                     IndexPermutation(o2pTrialDofs)));
    return result;
#else // without Ahmed
    throw std::runtime_error("To enable assembly in ACA mode, recompile BEM++ "
                             "with the symbol WITH_AHMED defined.");
#endif
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
AcaGlobalAssembler<ValueType>::assembleWeakForm(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        LocalAssembler& localAssembler,
        const AssemblyOptions& options)
{
    std::vector<LocalAssembler*> localAssemblers(1, &localAssembler);
    std::vector<const DiscreteLinOp*> sparseTermsToAdd;

    return assembleWeakForm(testSpace, trialSpace, localAssemblers,
                            sparseTermsToAdd, options);
}


#ifdef COMPILE_FOR_FLOAT
template class AcaGlobalAssembler<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class AcaGlobalAssembler<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class AcaGlobalAssembler<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class AcaGlobalAssembler<std::complex<double> >;
#endif

} // namespace Bempp
