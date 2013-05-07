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


#include "hypersingular_integral_operator.hpp"

#include "aca_global_assembler.hpp"
#include "assembly_options.hpp"
#include "dense_global_assembler.hpp"
#include "discrete_boundary_operator.hpp"
#include "context.hpp"
#include "local_assembler_construction_helper.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_operators.hpp"
#include "../fiber/quadrature_strategy.hpp"

#include "../common/boost_make_shared_fwd.hpp"

#include <stdexcept>
#include <iostream>

#include <tbb/tick_count.h>

namespace Bempp
{

template <typename BasisFunctionType, typename KernelType, typename ResultType>
HypersingularIntegralOperator<BasisFunctionType, KernelType, ResultType>::
HypersingularIntegralOperator(
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label,
        int symmetry) :
    Base(domain, range, dualToRange, label, symmetry)
{
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
bool
HypersingularIntegralOperator<BasisFunctionType, KernelType, ResultType>::isLocal() const
{
    return false;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::pair<
shared_ptr<typename HypersingularIntegralOperator<
BasisFunctionType, KernelType, ResultType>::LocalAssembler>,
shared_ptr<typename HypersingularIntegralOperator<
BasisFunctionType, KernelType, ResultType>::LocalAssembler>
>
HypersingularIntegralOperator<BasisFunctionType, KernelType, ResultType>::makeAssemblers(
        const QuadratureStrategy& quadStrategy,
        const AssemblyOptions& options) const
{
    typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
    typedef std::vector<const Fiber::Basis<BasisFunctionType>*> BasisPtrVector;

    const bool verbose = (options.verbosityLevel() >= VerbosityLevel::DEFAULT);

    shared_ptr<RawGridGeometry> testRawGeometry, trialRawGeometry;
    shared_ptr<GeometryFactory> testGeometryFactory, trialGeometryFactory;
    shared_ptr<Fiber::OpenClHandler> openClHandler;
    shared_ptr<BasisPtrVector> testBases, trialBases;
    bool cacheSingularIntegrals;

    if (verbose)
        std::cout << "Collecting data for assembler construction..." << std::endl;
       this->collectDataForAssemblerConstruction(options,
                                        testRawGeometry, trialRawGeometry,
                                        testGeometryFactory, trialGeometryFactory,
                                        testBases, trialBases,
                                        openClHandler, cacheSingularIntegrals);
    if (verbose)
        std::cout << "Data collection finished." << std::endl;

    bool makeSeparateOffDiagonalAssembler =
            options.assemblyMode() == AssemblyOptions::ACA &&
            options.acaOptions().globalAssemblyBeforeCompression == false;
    std::cout << "make? " << makeSeparateOffDiagonalAssembler << "\n";

    return reallyMakeAssemblers(quadStrategy,
                                testGeometryFactory, trialGeometryFactory,
                                testRawGeometry, trialRawGeometry,
                                testBases, trialBases, openClHandler,
                                options.parallelizationOptions(),
                                options.verbosityLevel(),
                                cacheSingularIntegrals,
                                makeSeparateOffDiagonalAssembler);
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::pair<
shared_ptr<typename HypersingularIntegralOperator<
BasisFunctionType, KernelType, ResultType>::LocalAssembler>,
shared_ptr<typename HypersingularIntegralOperator<
BasisFunctionType, KernelType, ResultType>::LocalAssembler>
>
HypersingularIntegralOperator<BasisFunctionType, KernelType, ResultType>::
reallyMakeAssemblers(
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
        bool cacheSingularIntegrals,
        bool makeSeparateOffDiagonalAssembler) const
{
    std::pair<shared_ptr<LocalAssembler>, shared_ptr<LocalAssembler> > result;
    // first element: "standard" assembler
    // second element: assembler used for admissible (off-diagonal)
    // H-matrix blocks in "disassembled mode"
    result.first.reset(quadStrategy.makeAssemblerForIntegralOperators(
                           testGeometryFactory, trialGeometryFactory,
                           testRawGeometry, trialRawGeometry,
                           testBases, trialBases,
                           make_shared_from_ref(testTransformations()),
                           make_shared_from_ref(kernels()),
                           make_shared_from_ref(trialTransformations()),
                           make_shared_from_ref(integral()),
                           openClHandler, parallelizationOptions, verbosityLevel,
                           cacheSingularIntegrals).release());
    if (makeSeparateOffDiagonalAssembler)
        result.second.reset(quadStrategy.makeAssemblerForIntegralOperators(
                                testGeometryFactory, trialGeometryFactory,
                                testRawGeometry, trialRawGeometry,
                                testBases, trialBases,
                                make_shared_from_ref(offDiagonalTestTransformations()),
                                make_shared_from_ref(offDiagonalKernels()),
                                make_shared_from_ref(offDiagonalTrialTransformations()),
                                make_shared_from_ref(offDiagonalIntegral()),
                                openClHandler, parallelizationOptions, verbosityLevel,
                                false /*cacheSingularIntegrals*/).release());
    else
        result.second = result.first;
    return result;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
HypersingularIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakFormImpl(
    const Context<BasisFunctionType, ResultType>& context) const
{
    bool verbose = (context.assemblyOptions().verbosityLevel() >=
                    VerbosityLevel::DEFAULT);
    if (verbose)
        std::cout << "Assembling the weak form of operator '"
                  << this->label() << "'..." << std::endl;

    tbb::tick_count start = tbb::tick_count::now();
    std::pair<shared_ptr<LocalAssembler>, shared_ptr<LocalAssembler> > assemblers =
        makeAssemblers(*context.quadStrategy(), context.assemblyOptions());
    shared_ptr<DiscreteBoundaryOperator<ResultType> > result =
        assembleWeakFormInternal(*assemblers.first, *assemblers.second, context);
    tbb::tick_count end = tbb::tick_count::now();

    if (verbose)
        std::cout << "Assembly of the weak form of operator '" << this->label()
                  << "' took " << (end - start).seconds() << " s" << std::endl;
    return result;
}

template <typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
HypersingularIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakFormInternal(
        LocalAssembler& standardAssembler, LocalAssembler& offDiagonalAssembler,
        const Context<BasisFunctionType, ResultType>& context) const
{
    switch (context.assemblyOptions().assemblyMode()) {
    case AssemblyOptions::DENSE:
        return shared_ptr<DiscreteBoundaryOperator<ResultType> >(
                    assembleWeakFormInDenseMode(standardAssembler, context).release());
    case AssemblyOptions::ACA:
        return shared_ptr<DiscreteBoundaryOperator<ResultType> >(
                    assembleWeakFormInAcaMode(
                        standardAssembler, offDiagonalAssembler, context).release());
    default:
        throw std::runtime_error(
                    "HypersingularIntegralOperator::assembleWeakFormInternalImpl(): "
                    "invalid assembly mode");
    }
}

// UNDOCUMENTED PRIVATE METHODS

/** \cond PRIVATE */

template <typename BasisFunctionType, typename KernelType, typename ResultType>
std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
HypersingularIntegralOperator<BasisFunctionType, KernelType, ResultType>::
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
HypersingularIntegralOperator<BasisFunctionType, KernelType, ResultType>::
assembleWeakFormInAcaMode(
        LocalAssembler& standardAssembler, LocalAssembler& offDiagonalAssembler,
        const Context<BasisFunctionType, ResultType>& context) const
{
    const Space<BasisFunctionType>& testSpace = *this->dualToRange();
    const Space<BasisFunctionType>& trialSpace = *this->domain();

    return AcaGlobalAssembler<BasisFunctionType, ResultType>::assembleDetachedWeakForm(
                testSpace, trialSpace, standardAssembler, offDiagonalAssembler,
                context, this->symmetry() & SYMMETRIC);
}
/** \endcond */

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_KERNEL_AND_RESULT(HypersingularIntegralOperator);

} // namespace Bempp
