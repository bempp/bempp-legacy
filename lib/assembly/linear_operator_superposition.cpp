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

#include "linear_operator_superposition.hpp"

#include "aca_global_assembler.hpp"
#include "assembly_options.hpp"
#include "discrete_dense_linear_operator.hpp"
#include "discrete_linear_operator_superposition.hpp"
#include "elementary_linear_operator.hpp"

#include "../common/auto_timer.hpp"
#include "../fiber/basis.hpp"
#include "../fiber/local_assembler_for_operators.hpp"
#include "../fiber/opencl_handler.hpp"
#include "../fiber/raw_grid_geometry.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/grid.hpp"
#include "../space/space.hpp"

namespace Bempp
{


template <typename ValueType>
LinearOperatorSuperposition<ValueType>::LinearOperatorSuperposition(
        const LinearOperator<ValueType>& term1,
        const LinearOperator<ValueType>& term2)
    : LinearOperator<ValueType>(term1.getTestSpace(),term1.getTrialSpace())
{

    this->addLocalOperatorsMultipliers(term1.getLocalOperators(),term1.getMultipliers());
    this->addLocalOperatorsMultipliers(term2.getLocalOperators(),term2.getMultipliers());
}

template <typename ValueType>
LinearOperatorSuperposition<ValueType>::LinearOperatorSuperposition(
        const LinearOperator<ValueType>& term,
        const ValueType& scalar)
: LinearOperator<ValueType>(term.getTestSpace(),term.getTrialSpace())
{
    const std::vector<ValueType>& m=term.getMultipliers();
    std::vector<ValueType> scaledMultipliers;
    for (int i=0;i<m.size();i++) scaledMultipliers.push_back(scalar*m[i]);
    this->addLocalOperatorsMultipliers(term.getLocalOperators(),scaledMultipliers);

}


template <typename ValueType>
int LinearOperatorSuperposition<ValueType>::trialComponentCount() const
{
    return this->getLocalOperators()[0]->trialComponentCount();
}

template <typename ValueType>
int LinearOperatorSuperposition<ValueType>::testComponentCount() const
{
    return this->getLocalOperators()[0]->testComponentCount();
}

template <typename ValueType>
bool LinearOperatorSuperposition<ValueType>::supportsRepresentation(
        AssemblyOptions::Representation repr) const
{
    return (repr == AssemblyOptions::DENSE || repr == AssemblyOptions::ACA);
}

template <typename ValueType>
std::auto_ptr<DiscreteLinearOperator<ValueType> >
LinearOperatorSuperposition<ValueType>::assembleWeakForm(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    switch (options.operatorRepresentation())
    {
    case AssemblyOptions::DENSE:
        return assembleWeakFormInDenseMode(factory, options);
    case AssemblyOptions::ACA:
        return assembleWeakFormInAcaMode(factory, options);
    default:
        return assembleWeakFormInArbitraryMode(factory, options);
    }
}

template <typename ValueType>
std::auto_ptr<DiscreteLinearOperator<ValueType> >
LinearOperatorSuperposition<ValueType>::assembleWeakFormInDenseMode(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    typedef DiscreteLinearOperator<ValueType> DiscreteLinOp;
    typedef DiscreteDenseLinearOperator<ValueType> DiscreteDenseLinOp;

    const std::vector<ElementaryLinearOperator<ValueType> const*> localOperators=this->getLocalOperators();
    const std::vector<ValueType> multipliers=this->getMultipliers();

    // Gather matrices of individual operators
    boost::ptr_vector<DiscreteLinOp> discreteOps;
    for (int i = 0; i < localOperators.size(); ++i)
    {
        std::auto_ptr<DiscreteLinOp> discreteOp =
                localOperators[i]->assembleWeakForm(factory, options);
        discreteOps.push_back(discreteOp);
    }

    // Add the matrices together
    arma::Mat<ValueType> sum;
    sum = discreteOps[0].asMatrix();
    for (int i = 1; i < discreteOps.size(); ++i)
            sum += discreteOps[i].asMatrix()*multipliers[i];

    return std::auto_ptr<DiscreteLinOp>(new DiscreteDenseLinOp(sum));
}

template <typename ValueType>
std::auto_ptr<DiscreteLinearOperator<ValueType> >
LinearOperatorSuperposition<ValueType>::assembleWeakFormInAcaMode(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    typedef DiscreteLinearOperator<ValueType> DiscreteLinOp;

    const Space<ValueType>& testSpace = this->getTestSpace();
    const Space<ValueType>& trialSpace = this->getTrialSpace();

    const std::vector<ElementaryLinearOperator<ValueType> const*> localOperators=this->getLocalOperators();
    const std::vector<ValueType> multipliers=this->getMultipliers();

    AutoTimer timer("\nAssembly took ");

    if (!testSpace.dofsAssigned() || !trialSpace.dofsAssigned())
        throw std::runtime_error(
                "LinearOperatorSuperposition::assembleWeakFormInAcaMode(): "
                "degrees of freedom must be assigned "
                "before calling assembleWeakForm()");
    if (&testSpace.grid() != &trialSpace.grid())
        throw std::runtime_error(
                "LinearOperatorSuperposition::assembleWeakFormInAcaMode(): "
                "testSpace and trialSpace must be defined over the same grid");

    // Prepare data for construction of local assembler

    const Grid& grid = trialSpace.grid();
    std::auto_ptr<GridView> view = grid.leafView();
    const int elementCount = view->entityCount(0);

    // REFACT The following two blocks might disappear in the constructor of
    // LocalAssemblerFactory

    // Gather geometric data
    Fiber::RawGridGeometry<ValueType> rawGeometry(grid.dim(), grid.dimWorld());
    view->getRawElementData(
                rawGeometry.vertices(), rawGeometry.elementCornerIndices(),
                rawGeometry.auxData());

    // Make geometry factory
    std::auto_ptr<GeometryFactory> geometryFactory =
            trialSpace.grid().elementGeometryFactory();

    // REFACT Basis retrieval might be moved into Space

    // Get pointers to test and trial bases of each element
    std::vector<const Fiber::Basis<ValueType>*> testBases;
    std::vector<const Fiber::Basis<ValueType>*> trialBases;
    testBases.reserve(elementCount);
    trialBases.reserve(elementCount);

    std::auto_ptr<EntityIterator<0> > it = view->entityIterator<0>();
    while (!it->finished())
    {
        const Entity<0>& element = it->entity();
        testBases.push_back(&testSpace.basis(element));
        trialBases.push_back(&trialSpace.basis(element));
        it->next();
    }

    // REFACT This will disappear in the constructor of LocalAssemblerFactory
    Fiber::OpenClHandler<ValueType, int> openClHandler(options.openClOptions());

    // REFACT This is unfortunately going to stay
    bool cacheSingularIntegrals =
            (options.singularIntegralCaching() == AssemblyOptions::YES ||
             (options.singularIntegralCaching() == AssemblyOptions::AUTO &&
              options.parallelism() == AssemblyOptions::OPEN_CL));


    // Construct local assemblers. Immediately assemble sparse terms in sparse
    // mode. Populate a vector of dense terms for subsequent ACA-mode assembly.
    boost::ptr_vector<DiscreteLinOp> sparseDiscreteTerms;
    boost::ptr_vector<LocalAssembler> denseTermLocalAssemblers;

    std::vector<ValueType> sparseTermsMultipliers;
    std::vector<ValueType> denseTermsMultipliers;

    for (int i = 0; i < localOperators.size(); ++i)
    {
        ElementaryLinearOperator<ValueType> const* term = localOperators[i];

        // Create local assembler for the current term
        std::auto_ptr<LocalAssembler> assembler = term->makeAssembler(
                    factory,
                    *geometryFactory, rawGeometry,
                    testBases, trialBases,
                    openClHandler, cacheSingularIntegrals);

        if (term->supportsRepresentation(AssemblyOptions::SPARSE))
        {
            std::auto_ptr<DiscreteLinOp> discreteTerm =
                    term->assembleWeakFormInternal(*assembler, options);
            sparseDiscreteTerms.push_back(discreteTerm);
            sparseTermsMultipliers.push_back(multipliers[i]);
        }
        else
        {
            denseTermLocalAssemblers.push_back(assembler);
            denseTermsMultipliers.push_back(multipliers[i]);
            assert(!assembler.get());
        }
    }

    // Convert boost::ptr_vectors to std::vectors
    std::vector<LocalAssembler*> stlDenseTermLocalAssemblers(
                denseTermLocalAssemblers.size());
    for (int i = 0; i < denseTermLocalAssemblers.size(); ++i)
        stlDenseTermLocalAssemblers[i] = &denseTermLocalAssemblers[i];

    std::vector<const DiscreteLinOp*> stlSparseDiscreteTerms(
                sparseDiscreteTerms.size());
    for (int i = 0; i < sparseDiscreteTerms.size(); ++i)
        stlSparseDiscreteTerms[i] = &sparseDiscreteTerms[i];

    // Assemble dense terms in ACA mode, simultaneously adding the sparse terms
    return AcaGlobalAssembler<ValueType>::assembleWeakForm(
                testSpace, trialSpace,
                stlDenseTermLocalAssemblers,
                stlSparseDiscreteTerms,
                denseTermsMultipliers,
                sparseTermsMultipliers,
                options);
}

template <typename ValueType>
std::auto_ptr<DiscreteLinearOperator<ValueType> >
LinearOperatorSuperposition<ValueType>::assembleWeakFormInArbitraryMode(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    // General (less efficient) implementation

    typedef DiscreteLinearOperator<ValueType> DiscreteLinOp;
    typedef DiscreteLinearOperatorSuperposition<ValueType>
            DiscreteSuperposition;

    const std::vector<ElementaryLinearOperator<ValueType> const*> localOperators=this->getLocalOperators();
    const std::vector<ValueType> multipliers=this->getMultipliers();


    boost::ptr_vector<DiscreteLinOp> discreteOps;
    for (int i = 0; i < localOperators.size(); ++i)
    {
        std::auto_ptr<DiscreteLinOp> discreteOp =
                localOperators[i]->assembleWeakForm(factory, options);
        discreteOps.push_back(discreteOp);
    }
    return std::auto_ptr<DiscreteLinOp>(new DiscreteSuperposition(discreteOps,this->getMultipliers()));
}


#ifdef COMPILE_FOR_FLOAT
template class LinearOperatorSuperposition<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class LinearOperatorSuperposition<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class LinearOperatorSuperposition<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class LinearOperatorSuperposition<std::complex<double> >;
#endif

} // namespace Bempp
