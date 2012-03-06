#include "elementary_linear_operator.hpp"

#include "assembly_options.hpp"
#include "discrete_dense_scalar_valued_linear_operator.hpp"
#include "discrete_vector_valued_linear_operator.hpp"

#include "../common/multidimensional_arrays.hpp"
#include "../common/not_implemented_error.hpp"
#include "../fiber/integration_manager.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../space/space.hpp"

#include <armadillo>
#include <boost/ptr_container/ptr_vector.hpp>
#include <stdexcept>

#ifdef WITH_AHMED
#include "ahmed_aux.hpp"
#include "discrete_aca_scalar_valued_linear_operator.hpp"
#include "weak_form_aca_assembly_helper.hpp"
#endif

// The spaces should be given the grid in constructor!

namespace Bempp
{

template <typename ValueType>
std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
ElementaryLinearOperator<ValueType>::assembleOperator(
        const arma::Mat<ctype>& testPoints,
        const Space<ValueType>& trialSpace,
        const IntegrationManagerFactory& factory,
        const AssemblyOptions& options) const
{
    if (!trialSpace.dofsAssigned())
        throw std::runtime_error("ElementaryLinearOperator::assembleOperator(): "
                                 "degrees of freedom must be assigned "
                                 "before calling assembleOperator()");

    arma::Mat<ValueType> vertices;
    arma::Mat<int> elementCorners;
    arma::Mat<char> auxData;
    trialSpace.grid().leafView()->getRawElementData(vertices, elementCorners,
                                                    auxData);
    std::auto_ptr<GeometryFactory> geometryFactory =
            trialSpace.grid().elementGeometryFactory();

    std::auto_ptr<IntegrationManager> intMgr =
            makeIntegrationManager(factory, *geometryFactory,
                                   vertices, elementCorners, auxData);

    switch (options.mode)
    {
    case ASSEMBLY_MODE_DENSE:
        return assembleOperatorInDenseMode(
                    testPoints, trialSpace, *intMgr, options);
    case ASSEMBLY_MODE_ACA:
        return assembleOperatorInAcaMode(
                    testPoints, trialSpace, *intMgr, options);
    case ASSEMBLY_MODE_FMM:
        throw std::runtime_error("ElementaryLinearOperator::assembleOperator(): "
                                 "assembly mode FMM is not implemented yet");
    default:
        throw std::runtime_error("ElementaryLinearOperator::assembleOperator(): "
                                 "invalid assembly mode");
    }
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
ElementaryLinearOperator<ValueType>::assembleWeakForm(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        const IntegrationManagerFactory& factory,
        const AssemblyOptions& options) const
{
    if (!testSpace.dofsAssigned() || !trialSpace.dofsAssigned())
        throw std::runtime_error("ElementaryLinearOperator::assembleWeakForm(): "
                                 "degrees of freedom must be assigned "
                                 "before calling assembleOperator()");
    if (&testSpace.grid() != &trialSpace.grid())
        throw std::runtime_error("ElementaryLinearOperator::assembleWeakForm(): "
                                 "testSpace and trialSpace must be defined over "
                                 "the same grid");

    arma::Mat<ctype> vertices;
    arma::Mat<int> elementCorners;
    arma::Mat<char> auxData;
    testSpace.grid().leafView()->getRawElementData(vertices, elementCorners,
                                                   auxData);
    std::auto_ptr<GeometryFactory> geometryFactory =
            testSpace.grid().elementGeometryFactory();

    std::auto_ptr<IntegrationManager > intMgr =
            makeIntegrationManager(factory, *geometryFactory,
                                   vertices, elementCorners, auxData);

    switch (options.mode)
    {
    case ASSEMBLY_MODE_DENSE:
        return assembleWeakFormInDenseMode(
                    testSpace, trialSpace, *intMgr, options);
    case ASSEMBLY_MODE_ACA:
        return assembleWeakFormInAcaMode(
                    testSpace, trialSpace, *intMgr, options);
    case ASSEMBLY_MODE_FMM:
        throw std::runtime_error("ElementaryLinearOperator::assembleWeakForm(): "
                                 "assembly mode FMM is not implemented yet");
    default:
        throw std::runtime_error("ElementaryLinearOperator::assembleWeakForm(): "
                                 "invalid assembly mode");
    }
}

template <typename ValueType>
std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
ElementaryLinearOperator<ValueType>::assembleOperatorInDenseMode(
        const arma::Mat<ctype>& testPoints,
        const Space<ValueType>& trialSpace,
        typename ElementaryLinearOperator<ValueType>::IntegrationManager& intMgr,
        const AssemblyOptions& options) const
{
    throw NotImplementedError("ElementaryLinearOperator::"
                              "assembleOperatorInDenseMode(): "
                              "not implemented yet");
}

template <typename ValueType>
std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
ElementaryLinearOperator<ValueType>::assembleOperatorInAcaMode(
        const arma::Mat<ctype>& testPoints,
        const Space<ValueType>& trialSpace,
        typename ElementaryLinearOperator<ValueType>::IntegrationManager& intMgr,
        const AssemblyOptions& options) const
{
    throw NotImplementedError("ElementaryLinearOperator::"
                              "assembleOperatorInAcaMode(): "
                              "not implemented yet");
}


template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
ElementaryLinearOperator<ValueType>::assembleWeakFormInDenseMode(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        typename ElementaryLinearOperator<ValueType>::IntegrationManager& intMgr,
        const AssemblyOptions& options) const
{
    // Get the grid's leaf view so that we can iterate over elements
    std::auto_ptr<GridView> view = trialSpace.grid().leafView();

    // Create EntityPointers to all elements.
    // For now we assume they fit in memory...
    const int elementCount = view->entityCount(0);
    boost::ptr_vector<EntityPointer<0> > elementsOwner;
    elementsOwner.reserve(elementCount);

    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished())
    {
        elementsOwner.push_back(it->frozen());
        it->next();
    }
    std::vector<const EntityPointer<0>*> elements(elementCount);
    for (int i = 0; i < elementCount; ++i)
        elements[i] = &elementsOwner[i];

    // Create the operator's matrix
    arma::Mat<ValueType> result(testSpace.globalDofCount(),
                                trialSpace.globalDofCount());
    std::vector<arma::Mat<ValueType> > localResult;

    // Storage of global DOF indices corresponding to local DOFs on single
    // test and trial elements
    std::vector<GlobalDofIndex> trialGlobalDofs;
    std::vector<GlobalDofIndex> testGlobalDofs;

    // Loop over trial elements
    for (int trialIndex = 0; trialIndex < elementCount; ++trialIndex)
    {
        const EntityPointer<0>* trialElement = elements[trialIndex];

        // Evaluate integrals over pairs of the current trial element and
        // all the test elements
        evaluateLocalWeakForms(TEST_TRIAL, elements, *trialElement, ALL_DOFS,
                               testSpace, trialSpace, intMgr, localResult);

        // Retrieve global DOF indices corresponding to the local DOFs of the
        // current trial element
        // (note: later we might optimise by storing explicitly only first DOF indices
        // for the vertex, edge and bubble families)
        trialSpace.globalDofs(trialElement->entity(), trialGlobalDofs);

        // Loop over test indices
        for (int testIndex = 0; testIndex < elementCount; ++testIndex)
        {
            // Retrieve global DOF indices corresponding to the local DOFs of the
            // current test element (note: this information might be cached)
            testSpace.globalDofs(elements[testIndex]->entity(), testGlobalDofs);

            // Add the integrals to appropriate entries in the operator's matrix
            for (int trialDof = 0; trialDof < trialGlobalDofs.size(); ++trialDof)
                for (int testDof = 0; testDof < testGlobalDofs.size(); ++testDof)
                result(testGlobalDofs[testDof], trialGlobalDofs[trialDof]) +=
                        localResult[testIndex](testDof, trialDof);
        }
    }

    // Create and return a discrete operator represented by the matrix that
    // has just been calculated
    return std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >(
                new DiscreteDenseScalarValuedLinearOperator<ValueType>(result));
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
ElementaryLinearOperator<ValueType>::assembleWeakFormInAcaMode(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        typename ElementaryLinearOperator<ValueType>::IntegrationManager& intMgr,
        const AssemblyOptions& options) const
{
#ifdef WITH_AHMED
    typedef AhmedDofWrapper<ValueType> AhmedDofType;
    typedef DiscreteScalarValuedLinearOperator<ValueType> DiscreteLinOp;
    typedef DiscreteAcaScalarValuedLinearOperator<ValueType,
            AhmedDofType, AhmedDofType > DiscreteAcaLinOp;

    // Get the grid's leaf view so that we can iterate over elements
    std::auto_ptr<GridView> view = trialSpace.grid().leafView();

    // const int elementCount = view.entityCount(0);
    const int trialDofCount = trialSpace.globalDofCount();
    const int testDofCount = testSpace.globalDofCount();

#ifndef NDEBUG
    std::cout << "Generating cluster trees ... " << std::flush;
#endif
    // o2p: map of original indices to permuted indices
    // p2o: map of permuted indices to original indices
    arma::Col<unsigned int> o2pTestDofs(testDofCount);
    arma::Col<unsigned int> p2oTestDofs(testDofCount);
    arma::Col<unsigned int> o2pTrialDofs(trialDofCount);
    arma::Col<unsigned int> p2oTrialDofs(trialDofCount);
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
                ahmedTestDofCenters, o2pTestDofs.memptr(),
                0, testDofCount);
    testClusterTree.createClusterTree(
                options.acaMinimumBlockSize,
                o2pTestDofs.memptr(), p2oTestDofs.memptr());
    bemcluster<const AhmedDofType> trialClusterTree(
                ahmedTrialDofCenters, o2pTrialDofs.memptr(),
                0, trialDofCount);
    trialClusterTree.createClusterTree(
                options.acaMinimumBlockSize,
                o2pTrialDofs.memptr(), p2oTrialDofs.memptr());

#ifndef NDEBUG
    std::cout << "Test cluster count: " << testClusterTree.getncl()
              << "\nTrial cluster count: " << trialClusterTree.getncl()
              << std::endl;
    std::cout << "o2pTest:\n" << o2pTestDofs << std::endl;
    std::cout << "p2oTest:\n" << p2oTestDofs << std::endl;
#endif

    typedef bemblcluster<AhmedDofType, AhmedDofType> DoubleCluster;
    std::auto_ptr<DoubleCluster> doubleClusterTree(
                new DoubleCluster(0, 0, testDofCount, trialDofCount));
    unsigned int blockCount = 0;
    doubleClusterTree->subdivide(&testClusterTree, &trialClusterTree,
                                 options.acaEta * options.acaEta,
                                 blockCount);

#ifndef NDEBUG
    std::cout << "Double cluster count: " << blockCount << std::endl;
#endif

    std::auto_ptr<DiscreteLinOp> result;

    // OpenMP implementation also possible

    WeakFormAcaAssemblyHelper<ValueType>
            helper(*this, *view, testSpace, trialSpace,
                   p2oTestDofs, p2oTrialDofs, intMgr, options);

    mblock<ValueType>** blocks = 0;
    allocmbls(doubleClusterTree.get(), blocks);
    matgen_sqntl(helper, doubleClusterTree.get(), doubleClusterTree.get(),
                 options.acaRecompress, options.acaEps,
                 options.acaMaximumRank, blocks);
    result = std::auto_ptr<DiscreteLinOp>(
                new DiscreteAcaLinOp(testDofCount, trialDofCount,
                                     doubleClusterTree, blocks));
    return result;
#else // without Ahmed
    throw std::runtime_error("To enable assembly in ACA mode, recompile BEM++ "
                             "with the symbol WITH_AHMED defined.");
#endif
}

#ifdef COMPILE_FOR_FLOAT
template class ElementaryLinearOperator<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class ElementaryLinearOperator<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class ElementaryLinearOperator<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class ElementaryLinearOperator<std::complex<double> >;
#endif

} // namespace Bempp
