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

#include "abstract_boundary_operator_superposition_base.hpp"

#include "context.hpp"
#include "abstract_boundary_operator_sum.hpp"
#include "aca_global_assembler.hpp"
#include "discrete_boundary_operator_sum.hpp"
#include "discrete_dense_boundary_operator.hpp"
#include "discrete_null_boundary_operator.hpp"
#include "elementary_integral_operator_base.hpp"
#include "scaled_abstract_boundary_operator.hpp"
#include "scaled_discrete_boundary_operator.hpp"

#include "../common/boost_make_shared_fwd.hpp"
#include "../common/boost_ptr_vector_fwd.hpp"
#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"

#include <tbb/tick_count.h>

namespace Bempp
{

namespace
{

template <typename BasisFunctionType, typename ResultType>
void decomposeBoundaryOperatorRecursively(
        const BoundaryOperator<BasisFunctionType, ResultType>& op,
        ResultType weight,
        std::vector<BoundaryOperator<BasisFunctionType, ResultType> >& joinableOps,
        std::vector<BoundaryOperator<BasisFunctionType, ResultType> >& nonjoinableOps,
        std::vector<ResultType>& joinableOpWeights,
        std::vector<ResultType>& nonjoinableOpWeights)
{
    // It would be possible to replace this function, and its stack of dynamic
    // cast, with a virtual method AbstractBoundaryOperator::decompose(...).
    // But I'm not sure whether it would really be more elegant: this method
    // would anyway need to take a BoundaryOperator object (because we want to
    // store a BO, not an ABO, in the *Ops vectors). So its implementation
    // would need the rather awkward check assert(op.abstractOperator().get()
    // == this).

    typedef BoundaryOperator<BasisFunctionType, ResultType> Op;
    typedef AbstractBoundaryOperator<BasisFunctionType, ResultType> AOp;
    typedef AbstractBoundaryOperatorSum<BasisFunctionType, ResultType> AOpSum;
    typedef ElementaryIntegralOperatorBase<BasisFunctionType, ResultType>
            ElemIntegralOp;

    typedef ScaledAbstractBoundaryOperator<BasisFunctionType, ResultType> ScaledAOp;

    shared_ptr<const AOp> abstractOp = op.abstractOperator();
    if (shared_ptr<const AOpSum> concreteOp =
            boost::dynamic_pointer_cast<const AOpSum>(abstractOp)) {
        Op term;
        // process term 1
        term = concreteOp->term1();
        if (!op.context()->assemblyOptions().isJointAssemblyEnabled() ||
                !term.context()->assemblyOptions().isJointAssemblyEnabled()) {
            nonjoinableOps.push_back(term);
            nonjoinableOpWeights.push_back(weight);
        } else
            decomposeBoundaryOperatorRecursively(
                        term, weight, joinableOps, nonjoinableOps,
                        joinableOpWeights, nonjoinableOpWeights);
        // process term 2
        term = concreteOp->term2();
        if (!op.context()->assemblyOptions().isJointAssemblyEnabled() ||
                !term.context()->assemblyOptions().isJointAssemblyEnabled()) {
            nonjoinableOps.push_back(term);
            nonjoinableOpWeights.push_back(weight);
        } else
            decomposeBoundaryOperatorRecursively(
                        term, weight, joinableOps, nonjoinableOps,
                        joinableOpWeights, nonjoinableOpWeights);
    } else if (shared_ptr<const ScaledAOp> concreteOp =
               boost::dynamic_pointer_cast<const ScaledAOp>(abstractOp)) {
        ResultType multiplier = concreteOp->multiplier();
        Op multiplicand = concreteOp->multiplicand();

        if (!op.context()->assemblyOptions().isJointAssemblyEnabled() ||
                !multiplicand.context()->assemblyOptions().isJointAssemblyEnabled()) {
            nonjoinableOps.push_back(multiplicand);
            nonjoinableOpWeights.push_back(weight * multiplier);
        } else
            decomposeBoundaryOperatorRecursively(
                        multiplicand, weight * multiplier,
                        joinableOps, nonjoinableOps,
                        joinableOpWeights, nonjoinableOpWeights);
    } else if (shared_ptr<const ElemIntegralOp> concreteOp =
               boost::dynamic_pointer_cast<const ElemIntegralOp>(abstractOp)) {
        if (!op.context()->assemblyOptions().isJointAssemblyEnabled()) {
            nonjoinableOps.push_back(op);
            nonjoinableOpWeights.push_back(weight);
        } else {
            joinableOps.push_back(op);
            joinableOpWeights.push_back(weight);
        }
    } else {
        nonjoinableOps.push_back(op);
        nonjoinableOpWeights.push_back(weight);
    }
}

template <typename BasisFunctionType, typename ResultType>
void decomposeBoundaryOperator(
        const BoundaryOperator<BasisFunctionType, ResultType>& op,
        std::vector<BoundaryOperator<BasisFunctionType, ResultType> >& joinableOps,
        std::vector<BoundaryOperator<BasisFunctionType, ResultType> >& nonjoinableOps,
        std::vector<ResultType>& joinableOpWeights,
        std::vector<ResultType>& nonjoinableOpWeights)
{
    joinableOps.clear();
    nonjoinableOps.clear();
    joinableOpWeights.clear();
    nonjoinableOpWeights.clear();
    decomposeBoundaryOperatorRecursively(
                op, static_cast<ResultType>(1.),
                joinableOps, nonjoinableOps,
                joinableOpWeights, nonjoinableOpWeights);
}

} // namespace

template <typename BasisFunctionType_, typename ResultType_>
AbstractBoundaryOperatorSuperpositionBase<BasisFunctionType_, ResultType_>::
AbstractBoundaryOperatorSuperpositionBase(
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label,
        int symmetry) :
    Base(domain, range, dualToRange, label, symmetry)
{
}

template <typename BasisFunctionType_, typename ResultType_>
shared_ptr<DiscreteBoundaryOperator<ResultType_> >
AbstractBoundaryOperatorSuperpositionBase<BasisFunctionType_, ResultType_>::
assembleWeakFormImpl(
        const Context<BasisFunctionType, ResultType>& context) const
{
    typedef BoundaryOperator<BasisFunctionType, ResultType> Op;
    typedef DiscreteBoundaryOperator<ResultType> DiscreteOp;
    typedef DiscreteBoundaryOperatorSum<ResultType> DiscreteOpSum;
    typedef DiscreteNullBoundaryOperator<ResultType> DiscreteNullOp;
    typedef ScaledDiscreteBoundaryOperator<ResultType> ScaledDiscreteOp;

    const ResultType one = 1.;

    bool verbose = (context.assemblyOptions().verbosityLevel() >=
                    VerbosityLevel::DEFAULT);

    Op thisOp(make_shared_from_ref(context), make_shared_from_ref(*this));
    std::vector<Op> joinableOps, nonjoinableOps;
    std::vector<ResultType> joinableOpWeights, nonjoinableOpWeights;
    decomposeBoundaryOperator(thisOp, joinableOps, nonjoinableOps,
                              joinableOpWeights, nonjoinableOpWeights);
    assert(joinableOps.size() == joinableOpWeights.size());
    assert(nonjoinableOps.size() == nonjoinableOpWeights.size());

    // std::cout << "joinable ops:\n";
    // for (size_t i = 0; i < joinableOps.size(); ++i)
    //     std::cout << joinableOps[i].label() << std::endl;
    // std::cout << "Nonjoinable ops:\n";
    // for (size_t i = 0; i < nonjoinableOps.size(); ++i)
    //     std::cout << nonjoinableOps[i].label() << std::endl;

    // Pass the list of the joinable operators to the appropriate
    // assemblejointOperatorWeakForm*() function. This function will
    // assemble some (possibly all, possibly none) of these operators and
    // return the result, at the same time removing these operators nad their
    // weights from the lists.
    shared_ptr<DiscreteOp> result;
    if (context.assemblyOptions().assemblyMode() == AssemblyOptions::DENSE)
        result = assembleJointOperatorWeakFormInDenseMode(
                    joinableOps, joinableOpWeights, verbose);
    else if (context.assemblyOptions().assemblyMode() == AssemblyOptions::ACA)
        result = assembleJointOperatorWeakFormInAcaMode(
                    context, joinableOps, joinableOpWeights);
    else
        throw std::invalid_argument(
            "AbstractBoundaryOperatorSuperpositionBase::"
            "assembleWeakFormImpl(): unknown assembly mode");
    // unprocessed operators will be treated as nonjoinable
    nonjoinableOps.insert(nonjoinableOps.end(),
                            joinableOps.begin(), joinableOps.end());
    nonjoinableOpWeights.insert(nonjoinableOpWeights.end(),
                                  joinableOpWeights.begin(),
                                  joinableOpWeights.end());

    if (!result) {
        // There were no joinable operators. Process as many nonjoinable
        // operators as necessary to get a *non-const* discrete operator.
        assert(nonjoinableOps.size() >= 1);
        if (nonjoinableOps.size() == 1)
            return boost::make_shared<ScaledDiscreteOp>(
                        nonjoinableOpWeights[0],
                        nonjoinableOps[0].weakForm());
        else {
            if (nonjoinableOpWeights[0] == one) {
                if (nonjoinableOpWeights[1] == one)
                    result = boost::make_shared<DiscreteOpSum>(
                                nonjoinableOps[0].weakForm(),
                                nonjoinableOps[1].weakForm());
                else
                    result = boost::make_shared<DiscreteOpSum>(
                                nonjoinableOps[0].weakForm(),
                                boost::make_shared<ScaledDiscreteOp>(
                                    nonjoinableOpWeights[1],
                                    nonjoinableOps[1].weakForm()));
                nonjoinableOps.erase(nonjoinableOps.begin());
                nonjoinableOps.erase(nonjoinableOps.begin());
            } else {
                result = boost::make_shared<ScaledDiscreteOp>(
                            nonjoinableOpWeights[0],
                            nonjoinableOps[0].weakForm());
                nonjoinableOps.erase(nonjoinableOps.begin());
            }
        }
    }

    for (size_t i = 0; i < nonjoinableOps.size(); ++i)
        if (nonjoinableOpWeights[i] == one)
            result = boost::make_shared<DiscreteOpSum>(
                        result, nonjoinableOps[i].weakForm());
        else
            result = boost::make_shared<DiscreteOpSum>(
                        result, boost::make_shared<ScaledDiscreteOp>(
                            nonjoinableOpWeights[i],
                            nonjoinableOps[i].weakForm()));
    return result;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
AbstractBoundaryOperatorSuperpositionBase<BasisFunctionType, ResultType>::
assembleJointOperatorWeakFormInDenseMode(
        std::vector<BoundaryOperator<BasisFunctionType, ResultType> >& ops,
        std::vector<ResultType>& opWeights, bool verbose) const
{
    typedef DiscreteBoundaryOperator<ResultType> DiscreteOp;
    typedef DiscreteDenseBoundaryOperator<ResultType> DiscreteDenseOp;
    typedef DiscreteBoundaryOperatorSum<ResultType> DiscreteOpSum;
    typedef ScaledDiscreteBoundaryOperator<ResultType> ScaledDiscreteOp;

    size_t opCount = ops.size();
    assert(opWeights.size() == opCount);
    if (opCount == 0)
        return shared_ptr<DiscreteOp>();

    shared_ptr<DiscreteDenseOp> discreteDenseOpSum;
    shared_ptr<DiscreteOp> discreteNondenseOpSum;
    for (size_t i = 0; i < opCount; ++i) {
        // here we don't call ops[i].weakForm() because we don't want the
        // discrete weak forms of constituent operators to be cached
        // (would defeat the purpose: storing only *one* dense matrix)
        shared_ptr<DiscreteOp> discreteOp =
                ops[i].abstractOperator()->assembleWeakForm(*ops[i].context());
        if (shared_ptr<DiscreteDenseOp> discreteDenseOp =
                boost::dynamic_pointer_cast<DiscreteDenseOp>(discreteOp)) {
            if (discreteDenseOpSum)
                discreteDenseOpSum = boost::make_shared<DiscreteDenseOp>(
                    discreteDenseOpSum->asMatrix() +
                    opWeights[i] * discreteDenseOp->asMatrix());
            else { // this is the first dense operator encountered
                if (opWeights[i] == static_cast<ResultType>(1.))
                    discreteDenseOpSum = discreteDenseOp;
                else
                    discreteDenseOpSum = boost::make_shared<DiscreteDenseOp>(
                        opWeights[i] * discreteDenseOp->asMatrix());
            }
        } else {
            if (opWeights[i] != static_cast<ResultType>(1.))
                discreteOp = boost::make_shared<ScaledDiscreteOp>(
                    opWeights[i], discreteOp);
            if (discreteNondenseOpSum)
                discreteNondenseOpSum = boost::make_shared<DiscreteOpSum>(
                    discreteNondenseOpSum, discreteOp);
            else
                discreteNondenseOpSum = discreteOp;
        }
    }

    ops.clear();
    opWeights.clear();

    if (discreteDenseOpSum) {
        if (discreteNondenseOpSum)
            return boost::make_shared<DiscreteOpSum>(
                discreteDenseOpSum, discreteNondenseOpSum);
        else
            return discreteDenseOpSum;
    } else
        return discreteNondenseOpSum;
}

template <typename BasisFunctionType_, typename ResultType_>
shared_ptr<DiscreteBoundaryOperator<ResultType_> >
AbstractBoundaryOperatorSuperpositionBase<BasisFunctionType_, ResultType_>::
assembleJointOperatorWeakFormInAcaMode(
        const Context<BasisFunctionType, ResultType>& context,
        std::vector<BoundaryOperator<BasisFunctionType, ResultType> >& ops,
        std::vector<ResultType>& opWeights) const
{
    typedef BoundaryOperator<BasisFunctionType, ResultType> Op;
    typedef DiscreteBoundaryOperator<ResultType> DiscreteOp;
    typedef DiscreteBoundaryOperatorSum<ResultType> DiscreteOpSum;
    typedef ElementaryIntegralOperatorBase<BasisFunctionType, ResultType>
            ElemIntegralOp;
    size_t opCount = ops.size();
    assert(opWeights.size() == opCount);

    bool verbose = (context.assemblyOptions().verbosityLevel() >=
                    VerbosityLevel::DEFAULT);

    // Split operators into local and nonlocal
    std::vector<Op> localOps, nonlocalOps;
    std::vector<ResultType> localOpWeights, nonlocalOpWeights;
    for (size_t i = 0; i < opCount; ++i)
        if (ops[i].abstractOperator()->isLocal()) {
            localOps.push_back(ops[i]);
            localOpWeights.push_back(opWeights[i]);
        } else {
            nonlocalOps.push_back(ops[i]);
            nonlocalOpWeights.push_back(opWeights[i]);
        }

    shared_ptr<DiscreteOp> nonlocalPart;
    if (!nonlocalOps.empty()) {
        std::string label;
        if (verbose) {
            // Prepare label
            for (size_t i = 0; i < nonlocalOps.size(); ++i)
                if (nonlocalOpWeights[i] == static_cast<ResultType>(1.))
                    label += "(" + nonlocalOps[i].label() + ") + ";
                else
                    label += toString(nonlocalOpWeights[i]) + " * (" +
                            nonlocalOps[i].label() + ") + ";
            label = label.substr(0, label.size() - 3); // remove the trailing plus

            std::cout << "Assembling the weak form of operator '"
                      << label << "'..." << std::endl;
        }
        tbb::tick_count start = tbb::tick_count::now();

        // Collect data used in the construction of all assemblers
        typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
        typedef std::vector<const Fiber::Shapeset<BasisFunctionType>*> ShapesetPtrVector;

        shared_ptr<RawGridGeometry> testRawGeometry, trialRawGeometry;
        shared_ptr<GeometryFactory> testGeometryFactory, trialGeometryFactory;
        shared_ptr<ShapesetPtrVector> testShapesets, trialShapesets;

        if (verbose)
            std::cout << "Collecting data for assembler construction..." << std::endl;
        this->collectOptionsIndependentDataForAssemblerConstruction(
                    testRawGeometry, trialRawGeometry,
                    testGeometryFactory, trialGeometryFactory,
                    testShapesets, trialShapesets);
        if (verbose)
            std::cout << "Data collection finished." << std::endl;

        // Construct assemblers and determine overall symmetry
        boost::ptr_vector<LocalAssembler> assemblersForNonlocalTerms;
        int symmetry = 0xfffffff;

        for (size_t i = 0; i < nonlocalOps.size(); ++i) {
            // All operators in the 'ops' list are expected to be joinable, i.e.
            // their abstract operators should be instances of
            // ElementaryIntegralOperatorBase (because it is this interface
            // that defines the makeAssembler() function
            shared_ptr<const ElemIntegralOp> elemOp =
                    boost::dynamic_pointer_cast<const ElemIntegralOp>(
                        nonlocalOps[i].abstractOperator());
            assert(elemOp);
            const AssemblyOptions& options =
                    nonlocalOps[i].context()->assemblyOptions();
            shared_ptr<Fiber::OpenClHandler> openClHandler;
            bool cacheSingularIntegrals;
            this->collectOptionsDependentDataForAssemblerConstruction(
                        options, testRawGeometry, trialRawGeometry,
                        openClHandler, cacheSingularIntegrals);

            // Create assembler for the current term
            std::unique_ptr<LocalAssembler> assembler = elemOp->makeAssembler(
                        *nonlocalOps[i].context()->quadStrategy(),
                        testGeometryFactory, trialGeometryFactory,
                        testRawGeometry, trialRawGeometry,
                        testShapesets, trialShapesets,
                        openClHandler,
                        options.parallelizationOptions(),
                        options.verbosityLevel(),
                        cacheSingularIntegrals);
            assemblersForNonlocalTerms.push_back(assembler.release());
            symmetry &= elemOp->symmetry();
        }

        // Convert boost::ptr_vectors to std::vectors
        std::vector<LocalAssembler*> stlAssemblersForNonlocalTerms(
                    assemblersForNonlocalTerms.size());
        for (size_t i = 0; i < assemblersForNonlocalTerms.size(); ++i)
            stlAssemblersForNonlocalTerms[i] = &assemblersForNonlocalTerms[i];

        std::vector<const DiscreteOp*> stlSparseDiscreteTerms; // empty; unused
        std::vector<ResultType> sparseTermMultipliers; // empty; unused

        // Assemble nonlocal terms in ACA mode
        nonlocalPart.reset(
                    AcaGlobalAssembler<BasisFunctionType, ResultType>::
                    assembleDetachedWeakForm(
                        *this->dualToRange(), *this->domain(),
                        stlAssemblersForNonlocalTerms,
                        stlAssemblersForNonlocalTerms,
                        stlSparseDiscreteTerms,
                        nonlocalOpWeights,
                        sparseTermMultipliers,
                        context, // We're using here the
                                 // options (esp. ACA options)
                                 // of the superposition operator
                        false /* no symmetry, for the moment */).release());
        tbb::tick_count end = tbb::tick_count::now();
        if (verbose)
            std::cout << "Assembly of the weak form of operator '" << label
                      << "' took " << (end - start).seconds() << " s" << std::endl;
    }

    // For the moment, we don't do anything special with identity operators,
    // but we might consider adding their contributions to the H-matrix
    ops = localOps;
    opWeights = localOpWeights;

    return nonlocalPart;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
        AbstractBoundaryOperatorSuperpositionBase);

} // namespace Bempp
