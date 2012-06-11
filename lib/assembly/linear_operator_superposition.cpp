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
    addConstituentOperators(
                term1.constituentOperators(), term1.constituentOperatorWeights());
    addConstituentOperators(
                term2.constituentOperators(), term2.constituentOperatorWeights());
}

template <typename BasisFunctionType, typename ResultType>
LinearOperatorSuperposition<BasisFunctionType, ResultType>::
LinearOperatorSuperposition(
        const LinearOperator<BasisFunctionType, ResultType>& term,
        const ResultType& scalar) :
    LinearOperator<BasisFunctionType, ResultType>(term.testSpace(), term.trialSpace())
{
    const std::vector<ResultType>& weights = term.constituentOperatorWeights();
    std::vector<ResultType> scaledWeigths;
    for (size_t i = 0; i < weights.size(); i++)
        scaledWeigths.push_back(scalar * weights[i]);
    addConstituentOperators(term.constituentOperators(), scaledWeigths);
}

template <typename BasisFunctionType, typename ResultType>
int LinearOperatorSuperposition<BasisFunctionType, ResultType>::
trialComponentCount() const
{
    return m_constituentOperators[0]->trialComponentCount();
}

template <typename BasisFunctionType, typename ResultType>
int LinearOperatorSuperposition<BasisFunctionType, ResultType>::
testComponentCount() const
{
    return m_constituentOperators[0]->testComponentCount();
}

template <typename BasisFunctionType, typename ResultType>
std::vector<ElementaryLinearOperator<BasisFunctionType, ResultType> const*>
LinearOperatorSuperposition<BasisFunctionType, ResultType>::
constituentOperators() const
{
    return m_constituentOperators;
}

template <typename BasisFunctionType, typename ResultType>
std::vector<ResultType>
LinearOperatorSuperposition<BasisFunctionType, ResultType>::
constituentOperatorWeights() const
{
    return m_constituentOperatorWeights;
}

template <typename BasisFunctionType, typename ResultType>
bool LinearOperatorSuperposition<BasisFunctionType, ResultType>::
supportsRepresentation(
        AssemblyOptions::Representation repr) const
{
    return (repr == AssemblyOptions::DENSE || repr == AssemblyOptions::ACA);
}

template <typename BasisFunctionType, typename ResultType>
void LinearOperatorSuperposition<BasisFunctionType, ResultType>::
addConstituentOperators(
        const std::vector<ElementaryLinearOperator<
        BasisFunctionType, ResultType> const*>& operators,
        const std::vector<ResultType>& weights)
{
    if (operators.size() != weights.size())
        throw std::invalid_argument("LinearOperatorSuperposition::"
                                    "addConstituentOperators(): "
                                    "argument lengths do not match");
    m_constituentOperators.insert(m_constituentOperators.end(),
                                  operators.begin(), operators.end());
    m_constituentOperatorWeights.insert(m_constituentOperatorWeights.end(),
                                        weights.begin(), weights.end());
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
LinearOperatorSuperposition<BasisFunctionType, ResultType>::
assembleDetachedWeakFormImpl(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options,
        Symmetry symmetry) const
{
    switch (options.operatorRepresentation()) {
    case AssemblyOptions::DENSE:
        return assembleDetachedWeakFormInDenseMode(factory, options, symmetry);
    case AssemblyOptions::ACA:
        // return assembleDetachedWeakFormInAcaMode(factory, options, symmetry);
    default:
        return assembleDetachedWeakFormInArbitraryMode(factory, options, symmetry);
    }
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
LinearOperatorSuperposition<BasisFunctionType, ResultType>::
assembleDetachedWeakFormInDenseMode(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options,
        Symmetry symmetry) const
{
    typedef DiscreteLinearOperator<ResultType> DiscreteLinOp;
    typedef DiscreteDenseLinearOperator<ResultType> DiscreteDenseLinOp;

    const std::vector<ElementaryLinearOperator<BasisFunctionType, ResultType> const*>
            operators = this->constituentOperators();
    const std::vector<ResultType>& weights = this->constituentOperatorWeights();

    // Gather matrices of individual operators
    boost::ptr_vector<DiscreteLinOp> discreteOps;
    for (size_t i = 0; i < operators.size(); ++i) {
        std::auto_ptr<DiscreteLinOp> discreteOp =
                operators[i]->assembleDetachedWeakForm(factory, options, symmetry);
        discreteOps.push_back(discreteOp);
    }

    // Add the matrices together
    arma::Mat<ResultType> sum;
    sum = discreteOps[0].asMatrix() * weights[0];
    for (size_t i = 1; i < discreteOps.size(); ++i) {
        sum += discreteOps[i].asMatrix() * weights[i];
    }

    return std::auto_ptr<DiscreteLinOp>(new DiscreteDenseLinOp(sum));
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
LinearOperatorSuperposition<BasisFunctionType, ResultType>::
assembleDetachedWeakFormInAcaMode(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options,
        Symmetry symmetry) const
{
    AutoTimer timer("\nAssembly took ");

    typedef DiscreteLinearOperator<ResultType> DiscreteLinOp;

    const std::vector<ElementaryLinearOperator<BasisFunctionType, ResultType> const*>
            operators = this->constituentOperators();
    const std::vector<ResultType>& weights = this->constituentOperatorWeights();

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

    for (size_t i = 0; i < operators.size(); ++i) {
        ElementaryLinearOperator<BasisFunctionType, ResultType> const* term =
                operators[i];

        // Create local assembler for the current term
        std::auto_ptr<LocalAssembler> assembler = term->makeAssembler(
                    factory,
                    testGeometryFactory, trialGeometryFactory,
                    testRawGeometry, trialRawGeometry,
                    testBases, trialBases,
                    openClHandler,
                    options.parallelisationOptions(), cacheSingularIntegrals);

        if (term->supportsRepresentation(AssemblyOptions::SPARSE)) {
            std::auto_ptr<DiscreteLinOp> discreteTerm =
                    term->assembleDetachedWeakFormInternal(*assembler, options);
            sparseDiscreteTerms.push_back(discreteTerm);
            sparseTermsMultipliers.push_back(weights[i]);
        } else {
            denseTermLocalAssemblers.push_back(assembler);
            denseTermsMultipliers.push_back(weights[i]);
            assert(!assembler.get());
        }
    }

    // Convert boost::ptr_vectors to std::vectors
    std::vector<LocalAssembler*> stlDenseTermLocalAssemblers(
                denseTermLocalAssemblers.size());
    for (size_t i = 0; i < denseTermLocalAssemblers.size(); ++i)
        stlDenseTermLocalAssemblers[i] = &denseTermLocalAssemblers[i];

    std::vector<const DiscreteLinOp*> stlSparseDiscreteTerms(
                sparseDiscreteTerms.size());
    for (size_t i = 0; i < sparseDiscreteTerms.size(); ++i)
        stlSparseDiscreteTerms[i] = &sparseDiscreteTerms[i];

    // Assemble dense terms in ACA mode, simultaneously adding the sparse terms
    return AcaGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
                this->testSpace(), this->trialSpace(),
                stlDenseTermLocalAssemblers,
                stlSparseDiscreteTerms,
                denseTermsMultipliers,
                sparseTermsMultipliers,
                options,
                symmetry);
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
LinearOperatorSuperposition<BasisFunctionType, ResultType>::
assembleDetachedWeakFormInArbitraryMode(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options,
        Symmetry symmetry) const
{
    // General (less efficient) implementation

    typedef DiscreteLinearOperator<ResultType> DiscreteLinOp;
    typedef DiscreteLinearOperatorSuperposition<ResultType>
            DiscreteSuperposition;

    const std::vector<ElementaryLinearOperator<BasisFunctionType, ResultType> const*>
            operators = this->constituentOperators();

    boost::ptr_vector<DiscreteLinOp> discreteOps;
    for (size_t i = 0; i < operators.size(); ++i) {
        std::auto_ptr<DiscreteLinOp> discreteOp =
                operators[i]->assembleDetachedWeakForm(factory, options, symmetry);
        discreteOps.push_back(discreteOp);
    }
    return std::auto_ptr<DiscreteLinOp>(
                new DiscreteSuperposition(discreteOps,
                                          this->constituentOperatorWeights()));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(LinearOperatorSuperposition);

} // namespace Bempp
