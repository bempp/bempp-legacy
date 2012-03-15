#include "elementary_integral_operator.hpp"

#include "assembly_options.hpp"
#include "discrete_dense_scalar_valued_linear_operator.hpp"
#include "discrete_vector_valued_linear_operator.hpp"

#include "../common/multidimensional_arrays.hpp"
#include "../common/not_implemented_error.hpp"
#include "../common/auto_timer.hpp"
#include "../fiber/local_assembler_factory.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../fiber/opencl_handler.hpp"
#include "../fiber/raw_grid_geometry.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../space/space.hpp"

#include <armadillo>
#include <boost/ptr_container/ptr_vector.hpp>
#include <stdexcept>

#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>
#include <tbb/task_scheduler_init.h>

#ifdef WITH_AHMED
#include "ahmed_aux.hpp"
#include "discrete_aca_scalar_valued_linear_operator.hpp"
#include "weak_form_aca_assembly_helper.hpp"
#endif

namespace Bempp
{

// Bodies of parallel loops

namespace
{

template <typename ValueType>
class DenseWeakFormAssemblerLoopBody
{
public:
    typedef tbb::spin_mutex MutexType;

    DenseWeakFormAssemblerLoopBody(const std::vector<int>& testIndices,
             const std::vector<std::vector<GlobalDofIndex> >& testGlobalDofs,
             const std::vector<std::vector<GlobalDofIndex> >& trialGlobalDofs,
             typename ElementaryIntegralOperator<ValueType>::LocalAssembler& assembler,
             arma::Mat<ValueType>& result, MutexType& mutex) :
        m_testIndices(testIndices),
        m_testGlobalDofs(testGlobalDofs), m_trialGlobalDofs(trialGlobalDofs),
        m_assembler(assembler), m_result(result), m_mutex(mutex)
    {}

    void operator() (const tbb::blocked_range<size_t>& r) const {
        const int elementCount = m_testIndices.size();
        std::vector<arma::Mat<ValueType> > localResult;
        for (size_t trialIndex = r.begin(); trialIndex != r.end(); ++trialIndex)
        {
            // Evaluate integrals over pairs of the current trial element and
            // all the test elements
            m_assembler.evaluateLocalWeakForms(TEST_TRIAL, m_testIndices, trialIndex,
                                               ALL_DOFS, localResult);

            const int trialDofCount = m_trialGlobalDofs[trialIndex].size();
            // Global assembly
            {
                MutexType::scoped_lock lock(m_mutex);
                // Loop over test indices
                for (int testIndex = 0; testIndex < elementCount; ++testIndex)
                {
                    const int testDofCount = m_testGlobalDofs[testIndex].size();
                    // Add the integrals to appropriate entries in the operator's matrix
                    for (int trialDof = 0; trialDof < trialDofCount; ++trialDof)
                        for (int testDof = 0; testDof < testDofCount; ++testDof)
                            m_result(m_testGlobalDofs[testIndex][testDof],
                                     m_trialGlobalDofs[trialIndex][trialDof]) +=
                                    localResult[testIndex](testDof, trialDof);
                }
            }
        }
    }

private:
    const std::vector<int>& m_testIndices;
    const std::vector<std::vector<GlobalDofIndex> >& m_trialGlobalDofs;
    const std::vector<std::vector<GlobalDofIndex> >& m_testGlobalDofs;
    // mutable OK because Assembler is thread-safe. (Alternative to "mutable" here:
    // make assembler's internal integrator map mutable)
    mutable typename ElementaryIntegralOperator<ValueType>::LocalAssembler& m_assembler;
    // mutable OK because write access to this matrix is protected by a mutex
    mutable arma::Mat<ValueType>& m_result;

    // mutex must be mutable because we need to lock and unlock it
    mutable MutexType& m_mutex;
};

} // namespace


template <typename ValueType>
std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
ElementaryIntegralOperator<ValueType>::assembleOperator(
        const arma::Mat<ctype>& testPoints,
        const Space<ValueType>& trialSpace,
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    if (!trialSpace.dofsAssigned())
        throw std::runtime_error("ElementaryIntegralOperator::assembleOperator(): "
                                 "degrees of freedom must be assigned "
                                 "before calling assembleOperator()");

    // Prepare local assembler

    std::auto_ptr<GridView> view = trialSpace.grid().leafView();
    const int elementCount = view->entityCount(0);

    // Gather geometric data
    Fiber::RawGridGeometry<ValueType> rawGeometry;
    view->getRawElementData(
                rawGeometry.vertices(), rawGeometry.elementCornerIndices(),
                rawGeometry.auxData());

    // Make geometry factory
    std::auto_ptr<GeometryFactory> geometryFactory =
            trialSpace.grid().elementGeometryFactory();

    // Get pointers to test and trial bases of each element
    std::vector<const Fiber::Basis<ValueType>*> trialBases;
    trialBases.reserve(elementCount);

    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished())
    {
        // const Entity<0>& element = it->entity();
        // TODO: make basis() accept const Entity<0>& instead of EntityPointer<0>
        trialBases.push_back(&trialSpace.basis(*it));
        it->next();
    }

    // Now create the assembler
    Fiber::OpenClHandler openClHandler(options.openClOptions());
    bool cacheSingularIntegrals =
            (options.singularIntegralCaching() == AssemblyOptions::YES ||
             (options.singularIntegralCaching() == AssemblyOptions::AUTO &&
              options.parallelism() == AssemblyOptions::OPEN_CL));

    std::auto_ptr<LocalAssembler> assembler =
            factory.make(*geometryFactory, rawGeometry,
                         trialBases,
                         kernel(), trialExpression(),
                         openClHandler, cacheSingularIntegrals);

    switch (options.operatorRepresentation())
    {
    case AssemblyOptions::DENSE:
        return assembleOperatorInDenseMode(
                    testPoints, trialSpace, *assembler, options);
    case AssemblyOptions::ACA:
        return assembleOperatorInAcaMode(
                    testPoints, trialSpace, *assembler, options);
    case AssemblyOptions::FMM:
        throw std::runtime_error("ElementaryIntegralOperator::assembleOperator(): "
                                 "assembly mode FMM is not implemented yet");
    default:
        throw std::runtime_error("ElementaryIntegralOperator::assembleOperator(): "
                                 "invalid assembly mode");
    }
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
ElementaryIntegralOperator<ValueType>::assembleWeakForm(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        const typename ElementaryIntegralOperator<ValueType>::LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    AutoTimer timer("\nAssembly took ");

    if (!testSpace.dofsAssigned() || !trialSpace.dofsAssigned())
        throw std::runtime_error("ElementaryIntegralOperator::assembleWeakForm(): "
                                 "degrees of freedom must be assigned "
                                 "before calling assembleOperator()");
    if (&testSpace.grid() != &trialSpace.grid())
        throw std::runtime_error("ElementaryIntegralOperator::assembleWeakForm(): "
                                 "testSpace and trialSpace must be defined over "
                                 "the same grid");

    // Prepare local assembler

    std::auto_ptr<GridView> view = trialSpace.grid().leafView();
    const int elementCount = view->entityCount(0);

    // Gather geometric data
    Fiber::RawGridGeometry<ValueType> rawGeometry;
    view->getRawElementData(
                rawGeometry.vertices(), rawGeometry.elementCornerIndices(),
                rawGeometry.auxData());

    // Make geometry factory
    std::auto_ptr<GeometryFactory> geometryFactory =
            trialSpace.grid().elementGeometryFactory();

    // Get pointers to test and trial bases of each element
    std::vector<const Fiber::Basis<ValueType>*> testBases;
    std::vector<const Fiber::Basis<ValueType>*> trialBases;
    testBases.reserve(elementCount);
    trialBases.reserve(elementCount);

    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished())
    {
        // const Entity<0>& element = it->entity();
        // TODO: make basis() accept const Entity<0>& instead of EntityPointer<0>
        testBases.push_back(&testSpace.basis(*it));
        trialBases.push_back(&trialSpace.basis(*it));
        it->next();
    }

    // Now create the assembler
    Fiber::OpenClHandler openClHandler(options.openClOptions());
    bool cacheSingularIntegrals =
            (options.singularIntegralCaching() == AssemblyOptions::YES ||
             (options.singularIntegralCaching() == AssemblyOptions::AUTO &&
              options.parallelism() == AssemblyOptions::OPEN_CL));

    std::auto_ptr<LocalAssembler> assembler =
            factory.make(*geometryFactory, rawGeometry,
                         testBases, trialBases,
                         testExpression(), kernel(), trialExpression(),
                         openClHandler, cacheSingularIntegrals);

    switch (options.operatorRepresentation())
    {
    case AssemblyOptions::DENSE:
        return assembleWeakFormInDenseMode(
                    testSpace, trialSpace, *assembler, options);
    case AssemblyOptions::ACA:
        return assembleWeakFormInAcaMode(
                    testSpace, trialSpace, *assembler, options);
    case AssemblyOptions::FMM:
        throw std::runtime_error("ElementaryIntegralOperator::assembleWeakForm(): "
                                 "assembly mode FMM is not implemented yet");
    default:
        throw std::runtime_error("ElementaryIntegralOperator::assembleWeakForm(): "
                                 "invalid assembly mode");
    }
}

template <typename ValueType>
std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
ElementaryIntegralOperator<ValueType>::assembleOperatorInDenseMode(
        const arma::Mat<ctype>& testPoints,
        const Space<ValueType>& trialSpace,
        typename ElementaryIntegralOperator<ValueType>::LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    throw NotImplementedError("ElementaryIntegralOperator::"
                              "assembleOperatorInDenseMode(): "
                              "not implemented yet");
}

template <typename ValueType>
std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
ElementaryIntegralOperator<ValueType>::assembleOperatorInAcaMode(
        const arma::Mat<ctype>& testPoints,
        const Space<ValueType>& trialSpace,
        typename ElementaryIntegralOperator<ValueType>::LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    throw NotImplementedError("ElementaryIntegralOperator::"
                              "assembleOperatorInAcaMode(): "
                              "not implemented yet");
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
ElementaryIntegralOperator<ValueType>::assembleWeakFormInDenseMode(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        typename ElementaryIntegralOperator<ValueType>::LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
    // Get the grid's leaf view so that we can iterate over elements
    std::auto_ptr<GridView> view = trialSpace.grid().leafView();
    const int elementCount = view->entityCount(0);

    // Global DOF indices corresponding to local DOFs on elements
    std::vector<std::vector<GlobalDofIndex> > trialGlobalDofs(elementCount);
    std::vector<std::vector<GlobalDofIndex> > testGlobalDofs(elementCount);

    // Gather global DOF lists
    const Mapper& mapper = view->elementMapper();
    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished())
    {
        const Entity<0>& element = it->entity();
        const int elementIndex = mapper.entityIndex(element);
        testSpace.globalDofs(element, testGlobalDofs[elementIndex]);
        trialSpace.globalDofs(element, trialGlobalDofs[elementIndex]);
        it->next();
    }

    // Make a vector of all element indices
    std::vector<int> testIndices(elementCount);
    for (int i = 0; i < elementCount; ++i)
        testIndices[i] = i;

    // Create the operator's matrix
    arma::Mat<ValueType> result(testSpace.globalDofCount(),
                                trialSpace.globalDofCount());
    result.fill(0.);

    typedef DenseWeakFormAssemblerLoopBody<ValueType> Body;
    typename Body::MutexType mutex;

    int maxThreadCount = 1;
    if (options.parallelism() == AssemblyOptions::TBB_AND_OPEN_MP)
    {
        if (options.maxThreadCount() == AssemblyOptions::AUTO)
            maxThreadCount = tbb::task_scheduler_init::automatic;
        else
            maxThreadCount = options.maxThreadCount();
    }
    tbb::task_scheduler_init scheduler(maxThreadCount);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, elementCount),
                      Body(testIndices, testGlobalDofs, trialGlobalDofs,
                           assembler, result, mutex));

//// Old serial code (TODO: decide whether to keep it behind e.g. #ifndef PARALLEL)
//    std::vector<arma::Mat<ValueType> > localResult;
//    // Loop over trial elements
//    for (int trialIndex = 0; trialIndex < elementCount; ++trialIndex)
//    {
//        // Evaluate integrals over pairs of the current trial element and
//        // all the test elements
//        assembler.evaluateLocalWeakForms(TEST_TRIAL, testIndices, trialIndex,
//                                         ALL_DOFS, localResult);

//        // Loop over test indices
//        for (int testIndex = 0; testIndex < elementCount; ++testIndex)
//            // Add the integrals to appropriate entries in the operator's matrix
//            for (int trialDof = 0; trialDof < trialGlobalDofs[trialIndex].size(); ++trialDof)
//                for (int testDof = 0; testDof < testGlobalDofs[testIndex].size(); ++testDof)
//                result(testGlobalDofs[testIndex][testDof],
//                       trialGlobalDofs[trialIndex][trialDof]) +=
//                        localResult[testIndex](testDof, trialDof);
//    }

    // Create and return a discrete operator represented by the matrix that
    // has just been calculated
    return std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >(
                new DiscreteDenseScalarValuedLinearOperator<ValueType>(result));
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
ElementaryIntegralOperator<ValueType>::assembleWeakFormInAcaMode(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        typename ElementaryIntegralOperator<ValueType>::LocalAssembler& assembler,
        const AssemblyOptions& options) const
{
#ifdef WITH_AHMED
    typedef AhmedDofWrapper<ValueType> AhmedDofType;
    typedef DiscreteScalarValuedLinearOperator<ValueType> DiscreteLinOp;
    typedef DiscreteAcaScalarValuedLinearOperator<ValueType,
            AhmedDofType, AhmedDofType > DiscreteAcaLinOp;

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

    // OpenMP implementation also possible

    WeakFormAcaAssemblyHelper<ValueType>
            helper(*this, *view, testSpace, trialSpace,
                   p2oTestDofs, p2oTrialDofs, assembler, options);

    mblock<ValueType>** blocks = 0;
    allocmbls(doubleClusterTree.get(), blocks);
    matgen_sqntl(helper, doubleClusterTree.get(), doubleClusterTree.get(),
                 acaOptions.recompress, acaOptions.eps,
                 acaOptions.maximumRank, blocks);
    result = std::auto_ptr<DiscreteLinOp>(
                new DiscreteAcaLinOp(testDofCount, trialDofCount,
                                     doubleClusterTree, blocks,
                                     o2pTestDofs, p2oTrialDofs));
    return result;
#else // without Ahmed
    throw std::runtime_error("To enable assembly in ACA mode, recompile BEM++ "
                             "with the symbol WITH_AHMED defined.");
#endif
}

#ifdef COMPILE_FOR_FLOAT
template class ElementaryIntegralOperator<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class ElementaryIntegralOperator<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class ElementaryIntegralOperator<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class ElementaryIntegralOperator<std::complex<double> >;
#endif

} // namespace Bempp
