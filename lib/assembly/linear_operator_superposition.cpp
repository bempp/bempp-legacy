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
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_operators.hpp"
#include "../fiber/opencl_handler.hpp"
#include "../fiber/raw_grid_geometry.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/grid.hpp"
#include "../space/space.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
LinearOperatorSuperposition<BasisFunctionType, ResultType>::
LinearOperatorSuperposition(
        const LinearOperator<BasisFunctionType, ResultType>& term1,
        const LinearOperator<BasisFunctionType, ResultType>& term2) :
    LinearOperator<BasisFunctionType, ResultType>(term1.testSpace(), term1.trialSpace())
{
    if (&term1.testSpace() != &term2.testSpace() ||
            &term1.trialSpace() != &term2.trialSpace())
        throw std::runtime_error(
                "LinearOperatorSuperposition::LinearOperatorSuperposition(): "
                "Spaces don't match");
    this->addLocalOperatorsAndMultipliers(
                term1.localOperators(), term1.multipliers());
    this->addLocalOperatorsAndMultipliers(
                term2.localOperators(), term2.multipliers());
}

template <typename BasisFunctionType, typename ResultType>
LinearOperatorSuperposition<BasisFunctionType, ResultType>::
LinearOperatorSuperposition(
        const LinearOperator<BasisFunctionType, ResultType>& term,
        const ResultType& scalar) :
    LinearOperator<BasisFunctionType, ResultType>(term.testSpace(), term.trialSpace())
{
    const std::vector<ResultType>& m = term.multipliers();
    std::vector<ResultType> scaledMultipliers;
    for (int i = 0; i < m.size(); i++)
        scaledMultipliers.push_back(scalar * m[i]);
    this->addLocalOperatorsAndMultipliers(
                term.localOperators(), scaledMultipliers);
}

template <typename BasisFunctionType, typename ResultType>
int LinearOperatorSuperposition<BasisFunctionType, ResultType>::
trialComponentCount() const
{
    return this->localOperators()[0]->trialComponentCount();
}

template <typename BasisFunctionType, typename ResultType>
int LinearOperatorSuperposition<BasisFunctionType, ResultType>::
testComponentCount() const
{
    return this->localOperators()[0]->testComponentCount();
}

template <typename BasisFunctionType, typename ResultType>
bool LinearOperatorSuperposition<BasisFunctionType, ResultType>::
supportsRepresentation(
        AssemblyOptions::Representation repr) const
{
    return (repr == AssemblyOptions::DENSE || repr == AssemblyOptions::ACA);
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
LinearOperatorSuperposition<BasisFunctionType, ResultType>::
assembleWeakForm(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    switch (options.operatorRepresentation()) {
    case AssemblyOptions::DENSE:
        return assembleWeakFormInDenseMode(factory, options);
    case AssemblyOptions::ACA:
        return assembleWeakFormInAcaMode(factory, options);
    default:
        return assembleWeakFormInArbitraryMode(factory, options);
    }
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
LinearOperatorSuperposition<BasisFunctionType, ResultType>::
assembleWeakFormInDenseMode(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    typedef DiscreteLinearOperator<ResultType> DiscreteLinOp;
    typedef DiscreteDenseLinearOperator<ResultType> DiscreteDenseLinOp;

    const std::vector<ElementaryLinearOperator<BasisFunctionType, ResultType> const*>
            localOperators = this->localOperators();
    const std::vector<ResultType>& multipliers = this->multipliers();

    // Gather matrices of individual operators
    boost::ptr_vector<DiscreteLinOp> discreteOps;
    for (int i = 0; i < localOperators.size(); ++i) {
        std::auto_ptr<DiscreteLinOp> discreteOp =
                localOperators[i]->assembleWeakForm(factory, options);
        discreteOps.push_back(discreteOp);
    }

    // Add the matrices together
    arma::Mat<ResultType> sum;
    sum = discreteOps[0].asMatrix() * multipliers[0];
    for (int i = 1; i < discreteOps.size(); ++i) {
        sum += discreteOps[i].asMatrix() * multipliers[i];
    }

    return std::auto_ptr<DiscreteLinOp>(new DiscreteDenseLinOp(sum));
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
LinearOperatorSuperposition<BasisFunctionType, ResultType>::
assembleWeakFormInAcaMode(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    AutoTimer timer("\nAssembly took ");

    typedef DiscreteLinearOperator<ResultType> DiscreteLinOp;

    const std::vector<ElementaryLinearOperator<BasisFunctionType, ResultType> const*>
            localOperators = this->localOperators();
    const std::vector<ResultType>& multipliers = this->multipliers();

    typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
    typedef std::vector<const Fiber::Basis<BasisFunctionType>*> BasisPtrVector;

    shared_ptr<RawGridGeometry> testRawGeometry, trialRawGeometry;
    shared_ptr<GeometryFactory> testGeometryFactory, trialGeometryFactory;
    shared_ptr<Fiber::OpenClHandler> openClHandler;
    shared_ptr<BasisPtrVector> testBases, trialBases;
    bool cacheSingularIntegrals;

    collectDataForAssemblerConstruction(options,
                                        testRawGeometry, trialRawGeometry,
                                        testGeometryFactory, trialGeometryFactory,
                                        testBases, trialBases,
                                        openClHandler, cacheSingularIntegrals);

    // Construct local assemblers. Immediately assemble sparse terms in sparse
    // mode. Populate a vector of dense terms for subsequent ACA-mode assembly.
    boost::ptr_vector<DiscreteLinOp> sparseDiscreteTerms;
    boost::ptr_vector<LocalAssembler> denseTermLocalAssemblers;

    std::vector<ResultType> sparseTermsMultipliers;
    std::vector<ResultType> denseTermsMultipliers;

    for (int i = 0; i < localOperators.size(); ++i) {
        ElementaryLinearOperator<BasisFunctionType, ResultType> const* term =
                localOperators[i];

        // Create local assembler for the current term
        std::auto_ptr<LocalAssembler> assembler = term->makeAssembler(
                    factory,
                    testGeometryFactory, testRawGeometry,
                    // TODO: add parameters for trial*
                    testBases, trialBases,
                    openClHandler,
                    options.parallelisationOptions(), cacheSingularIntegrals);

        if (term->supportsRepresentation(AssemblyOptions::SPARSE)) {
            std::auto_ptr<DiscreteLinOp> discreteTerm =
                    term->assembleWeakFormInternal(*assembler, options);
            sparseDiscreteTerms.push_back(discreteTerm);
            sparseTermsMultipliers.push_back(multipliers[i]);
        } else {
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
    return AcaGlobalAssembler<BasisFunctionType, ResultType>::assembleWeakForm(
                this->testSpace(), this->trialSpace(),
                stlDenseTermLocalAssemblers,
                stlSparseDiscreteTerms,
                denseTermsMultipliers,
                sparseTermsMultipliers,
                options);
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
LinearOperatorSuperposition<BasisFunctionType, ResultType>::assembleWeakFormInArbitraryMode(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    // General (less efficient) implementation

    typedef DiscreteLinearOperator<ResultType> DiscreteLinOp;
    typedef DiscreteLinearOperatorSuperposition<ResultType>
            DiscreteSuperposition;

    const std::vector<ElementaryLinearOperator<BasisFunctionType, ResultType> const*>
            localOperators = this->localOperators();

    boost::ptr_vector<DiscreteLinOp> discreteOps;
    for (int i = 0; i < localOperators.size(); ++i) {
        std::auto_ptr<DiscreteLinOp> discreteOp =
                localOperators[i]->assembleWeakForm(factory, options);
        discreteOps.push_back(discreteOp);
    }
    return std::auto_ptr<DiscreteLinOp>(
                new DiscreteSuperposition(discreteOps, this->multipliers()));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(LinearOperatorSuperposition);

} // namespace Bempp
