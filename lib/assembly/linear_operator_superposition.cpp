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
#include "discrete_dense_scalar_valued_linear_operator.hpp"
#include "discrete_scalar_valued_linear_operator_superposition.hpp"
#include "discrete_vector_valued_linear_operator_superposition.hpp"
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
        boost::ptr_vector<ElementaryLinearOperator<ValueType> >& terms)
{
   init(terms);
}

template <typename ValueType>
LinearOperatorSuperposition<ValueType>::LinearOperatorSuperposition(
        const boost::tuple<ElementaryLinearOperator<ValueType>*,
        ElementaryLinearOperator<ValueType>*>& terms)
{
    boost::ptr_vector<ElementaryLinearOperator<ValueType> > vTerms;
    vTerms.push_back(terms.template get<0>());
    vTerms.push_back(terms.template get<1>());
    init(vTerms);
}

template <typename ValueType>
int LinearOperatorSuperposition<ValueType>::trialComponentCount() const
{
    if (m_terms.empty())
        return 1;
    else
        return m_terms[0].trialComponentCount();
}

template <typename ValueType>
int LinearOperatorSuperposition<ValueType>::testComponentCount() const
{
    if (m_terms.empty())
        return 1;
    else
        return m_terms[0].testComponentCount();
}

template <typename ValueType>
bool LinearOperatorSuperposition<ValueType>::supportsRepresentation(
        AssemblyOptions::Representation repr) const
{
    return (repr == AssemblyOptions::DENSE || repr == AssemblyOptions::ACA);
}

template <typename ValueType>
void LinearOperatorSuperposition<ValueType>::init(
        boost::ptr_vector<ElementaryLinearOperator<ValueType> >& terms)
{
    if (!terms.empty())
    {
        int testComponentCount_ = terms[0].testComponentCount();
        int trialComponentCount_ = terms[0].trialComponentCount();
        for (int i = 0; i < terms.size(); ++i)
            if (testComponentCount_ != terms[0].testComponentCount() ||
                    trialComponentCount_ != terms[0].trialComponentCount())
                throw std::invalid_argument("LinearOperatorSuperposition::init(): "
                                            "incompatible operator dimensions");
        const int origTermCount = terms.size();
        m_terms.transfer(m_terms.end(), terms);
        const int newTermCount = m_terms.size();
        assert(origTermCount == newTermCount); // documentation of transfer() is
                                               // somewhat difficult to follow...

    }
}

template <typename ValueType>
std::auto_ptr<DiscreteVectorValuedLinearOperator<ValueType> >
LinearOperatorSuperposition<ValueType>::assembleOperator(
        const arma::Mat<ctype>& testPoints,
        const Space<ValueType>& trialSpace,
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    typedef DiscreteVectorValuedLinearOperator<ValueType> DiscreteLinOp;
    typedef DiscreteVectorValuedLinearOperatorSuperposition<ValueType>
            DiscreteSuperposition;

    boost::ptr_vector<DiscreteLinOp> discreteOps;
    for (int i = 0; i < m_terms.size(); ++i)
    {
        std::auto_ptr<DiscreteLinOp> discreteOp =
                m_terms[i].assembleOperator(testPoints, trialSpace, factory, options);
        discreteOps.push_back(discreteOp);
    }
    return std::auto_ptr<DiscreteLinOp>(new DiscreteSuperposition(discreteOps));
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
LinearOperatorSuperposition<ValueType>::assembleWeakForm(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    switch (options.operatorRepresentation())
    {
    case AssemblyOptions::DENSE:
        return assembleWeakFormInDenseMode(
                    testSpace, trialSpace, factory, options);
    case AssemblyOptions::ACA:
        return assembleWeakFormInAcaMode(
                    testSpace, trialSpace, factory, options);
    default:
        return assembleWeakFormInArbitraryMode(
                    testSpace, trialSpace, factory, options);
    }
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
LinearOperatorSuperposition<ValueType>::assembleWeakFormInDenseMode(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    typedef DiscreteScalarValuedLinearOperator<ValueType> DiscreteLinOp;
    typedef DiscreteDenseScalarValuedLinearOperator<ValueType> DiscreteDenseLinOp;

    // Gather matrices of individual operators
    boost::ptr_vector<DiscreteLinOp> discreteOps;
    for (int i = 0; i < m_terms.size(); ++i)
    {
        std::auto_ptr<DiscreteLinOp> discreteOp =
                m_terms[i].assembleWeakForm(testSpace, trialSpace, factory, options);
        discreteOps.push_back(discreteOp);
    }

    // Add the matrices together
    arma::Mat<ValueType> sum;
    if (discreteOps.empty())
    {
        sum.set_size(testSpace.globalDofCount(), trialSpace.globalDofCount());
        sum.fill(0.);
    }
    else
    {
        sum = discreteOps[0].asMatrix();
        for (int i = 1; i < discreteOps.size(); ++i)
            sum += discreteOps[i].asMatrix();
    }

    return std::auto_ptr<DiscreteLinOp>(new DiscreteDenseLinOp(sum));
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
LinearOperatorSuperposition<ValueType>::assembleWeakFormInAcaMode(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    typedef DiscreteScalarValuedLinearOperator<ValueType> DiscreteLinOp;

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

    std::auto_ptr<GridView> view = trialSpace.grid().leafView();
    const int elementCount = view->entityCount(0);

    // REFACT The following two blocks might disappear in the constructor of
    // LocalAssemblerFactory

    // Gather geometric data
    Fiber::RawGridGeometry<ValueType> rawGeometry;
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

    AssemblyOptions sparseModeOptions = options;
    sparseModeOptions.switchToSparse();

    // Construct local assemblers. Immediately assemble sparse terms in sparse
    // mode. Populate a vector of dense terms for subsequent ACA-mode assembly.
    boost::ptr_vector<DiscreteLinOp> sparseDiscreteTerms;
    boost::ptr_vector<LocalAssembler> denseTermLocalAssemblers;
    for (int i = 0; i < m_terms.size(); ++i)
    {
        const ElementaryLinearOperator<ValueType>& term = m_terms[i];

        // Create local assembler for the current term
        std::auto_ptr<LocalAssembler> assembler = term.makeAssembler(
                    factory,
                    *geometryFactory, rawGeometry,
                    testBases, trialBases,
                    openClHandler, cacheSingularIntegrals);

        if (term.supportsRepresentation(AssemblyOptions::SPARSE))
        {
            std::auto_ptr<DiscreteLinOp> discreteTerm =
                    term.assembleWeakFormInternal(
                        testSpace, trialSpace, *assembler, sparseModeOptions);
            sparseDiscreteTerms.push_back(discreteTerm);
        }
        else
        {
            denseTermLocalAssemblers.push_back(assembler);
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
                stlSparseDiscreteTerms, options);
}

template <typename ValueType>
std::auto_ptr<DiscreteScalarValuedLinearOperator<ValueType> >
LinearOperatorSuperposition<ValueType>::assembleWeakFormInArbitraryMode(
        const Space<ValueType>& testSpace,
        const Space<ValueType>& trialSpace,
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options) const
{
    // General (less efficient) implementation

    typedef DiscreteScalarValuedLinearOperator<ValueType> DiscreteLinOp;
    typedef DiscreteScalarValuedLinearOperatorSuperposition<ValueType>
            DiscreteSuperposition;

    boost::ptr_vector<DiscreteLinOp> discreteOps;
    for (int i = 0; i < m_terms.size(); ++i)
    {
        std::auto_ptr<DiscreteLinOp> discreteOp =
                m_terms[i].assembleWeakForm(testSpace, trialSpace, factory, options);
        discreteOps.push_back(discreteOp);
    }
    return std::auto_ptr<DiscreteLinOp>(new DiscreteSuperposition(discreteOps));
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
