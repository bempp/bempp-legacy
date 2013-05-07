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


#include "elementary_integral_operator.hpp"

#include "aca_global_assembler.hpp"
#include "assembly_options.hpp"
#include "dense_global_assembler.hpp"
#include "discrete_boundary_operator.hpp"
#include "context.hpp"
#include "local_assembler_construction_helper.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/quadrature_strategy.hpp"

#include "../common/boost_make_shared_fwd.hpp"

#include <stdexcept>
#include <iostream>

#include <tbb/tick_count.h>

namespace Bempp
{

template <typename BasisFunctionType, typename KernelType, typename ResultType>
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
ElementaryIntegralOperator(const shared_ptr<const Space<BasisFunctionType> >& domain,
                           const shared_ptr<const Space<BasisFunctionType> >& range,
                           const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
                           const std::string& label,
                           int symmetry) :
    Base(domain, range, dualToRange, label, symmetry)
{
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
bool
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::isLocal() const
{
    return false;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<typename ElementaryIntegralOperator<
BasisFunctionType, KernelType, ResultType>::LocalAssembler>
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
makeAssemblerImpl(
        const QuadratureStrategy& quadStrategy,
        const shared_ptr<const GeometryFactory>& testGeometryFactory,
        const shared_ptr<const GeometryFactory>& trialGeometryFactory,
        const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
        const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
        const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
        const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
        const shared_ptr<const Fiber::OpenClHandler>& openClHandler,
        const ParallelizationOptions& parallelizationOptions,
        VerbosityLevel::Level verbosityLevel,
        bool cacheSingularIntegrals) const
{
    return quadStrategy.makeAssemblerForIntegralOperators(
                testGeometryFactory, trialGeometryFactory,
                testRawGeometry, trialRawGeometry,
                testBases, trialBases,
                make_shared_from_ref(testTransformations()),
                make_shared_from_ref(kernels()),
                make_shared_from_ref(trialTransformations()),
                make_shared_from_ref(integral()),
                openClHandler, parallelizationOptions, verbosityLevel,
                cacheSingularIntegrals);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakFormImpl(
    const Context<BasisFunctionType, ResultType>& context) const
{
    bool verbose = (context.assemblyOptions().verbosityLevel() >=
                    VerbosityLevel::DEFAULT);
    if (verbose)
        std::cout << "Assembling the weak form of operator '"
                  << this->label() << "'..." << std::endl;

    tbb::tick_count start = tbb::tick_count::now();
    std::auto_ptr<LocalAssembler> assembler =
        this->makeAssembler(*context.quadStrategy(), context.assemblyOptions());
    shared_ptr<DiscreteBoundaryOperator<ResultType> > result =
        assembleWeakFormInternalImpl2(*assembler, context);
    tbb::tick_count end = tbb::tick_count::now();

    if (verbose)
        std::cout << "Assembly of the weak form of operator '" << this->label()
                  << "' took " << (end - start).seconds() << " s" << std::endl;
    return result;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakFormInternalImpl2(
        LocalAssembler& assembler,
        const Context<BasisFunctionType, ResultType>& context) const
{
    switch (context.assemblyOptions().assemblyMode()) {
    case AssemblyOptions::DENSE:
        return shared_ptr<DiscreteBoundaryOperator<ResultType> >(
                    assembleWeakFormInDenseMode(assembler, context).release());
    case AssemblyOptions::ACA:
        return shared_ptr<DiscreteBoundaryOperator<ResultType> >(
                    assembleWeakFormInAcaMode(assembler, context).release());
    default:
        throw std::runtime_error(
                    "ElementaryIntegralOperator::assembleWeakFormInternalImpl(): "
                    "invalid assembly mode");
    }
}

// UNDOCUMENTED PRIVATE METHODS

/** \cond PRIVATE */

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakFormInDenseMode(
        LocalAssembler& assembler,
        const Context<BasisFunctionType, ResultType>& context) const
{
    const Space<BasisFunctionType>& testSpace = *this->dualToRange();
    const Space<BasisFunctionType>& trialSpace = *this->domain();

    return DenseGlobalAssembler<BasisFunctionType, ResultType>::
            assembleDetachedWeakForm(testSpace, trialSpace, assembler, context);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakFormInAcaMode(
        LocalAssembler& assembler,
        const Context<BasisFunctionType, ResultType>& context) const
{
    const Space<BasisFunctionType>& testSpace = *this->dualToRange();
    const Space<BasisFunctionType>& trialSpace = *this->domain();

    // TODO: replace second assembler with assembler for admissible blocks
    return AcaGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
                testSpace, trialSpace, assembler, assembler,
                context, this->symmetry() & SYMMETRIC);
}
/** \endcond */

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(ElementaryIntegralOperator);

} // namespace Bempp
