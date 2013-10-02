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

#include "modified_helmholtz_3d_synthetic_boundary_operator_builder.hpp"

#include "abstract_boundary_operator.hpp"
#include "boundary_operator.hpp"
#include "context.hpp"
#include "synthetic_nonhypersingular_integral_operator_builder.hpp"

#include "../fiber/explicit_instantiation.hpp"

#include "../space/space.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename KernelType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
modifiedHelmholtz3dSyntheticBoundaryOperator(
    BoundaryOperator<BasisFunctionType, ResultType>
    (*constructor)(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& /*context*/,
        const shared_ptr<const Space<BasisFunctionType> >& /*domain*/,
        const shared_ptr<const Space<BasisFunctionType> >& /*range*/,
        const shared_ptr<const Space<BasisFunctionType> >& /*dualToRange*/,
        KernelType /*waveNumber*/,
        const std::string& /*label*/,
        int /*symmetry*/,
        bool /*useInterpolation*/,
        int /*interpPtsPerWavelength*/),
    const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
    const shared_ptr<const Space<BasisFunctionType> >& domain,
    const shared_ptr<const Space<BasisFunctionType> >& range,
    const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
    KernelType waveNumber,
    std::string label,
    int internalSymmetry,
    bool useInterpolation,
    int interpPtsPerWavelength,
    int maximumSyntheseSymmetry)
{
    if (!domain || !range || !dualToRange)
        throw std::invalid_argument(
            "modifiedHelmholtz3dSyntheticBoundaryOperator(): "
            "domain, range and dualToRange must not be null");

    shared_ptr<const Space<BasisFunctionType> > newDomain = domain;
    shared_ptr<const Space<BasisFunctionType> > newDualToRange = dualToRange;

    bool isBarycentric = (domain->isBarycentric() || dualToRange->isBarycentric());

    if (isBarycentric) {
        newDomain = domain->barycentricSpace(domain);
        newDualToRange = dualToRange->barycentricSpace(dualToRange);
    }

    AssemblyOptions internalAssemblyOptions = context->assemblyOptions();
    AcaOptions internalAcaOptions = internalAssemblyOptions.acaOptions();
    internalAcaOptions.mode = AcaOptions::GLOBAL_ASSEMBLY;
    internalAssemblyOptions.switchToAcaMode(internalAcaOptions);
    typedef Context<BasisFunctionType, ResultType> Ctx;
    shared_ptr<Ctx> internalContext(new Ctx(
            context->quadStrategy(), internalAssemblyOptions));
    shared_ptr<const Space<BasisFunctionType> > internalTrialSpace =
        newDomain->discontinuousSpace(domain);
    shared_ptr<const Space<BasisFunctionType> > internalTestSpace =
        newDualToRange->discontinuousSpace(dualToRange);
    if (label.empty())
        label = AbstractBoundaryOperator<BasisFunctionType, ResultType>::
            uniqueLabel();
    int syntheseSymmetry =
        (newDomain == newDualToRange && internalTrialSpace == internalTestSpace) ?
        maximumSyntheseSymmetry : 0;
    BoundaryOperator<BasisFunctionType, ResultType> internalOp =
        constructor(
            internalContext, internalTrialSpace, range /* or whatever */,
            internalTestSpace,
            waveNumber, "(" + label + ")_internal", internalSymmetry,
            useInterpolation, interpPtsPerWavelength);
    return syntheticNonhypersingularIntegralOperator(
            internalOp, newDomain, range, newDualToRange,
            internalTrialSpace, internalTestSpace,
            label, syntheseSymmetry);
}

#define INSTANTIATE_FUNCTION(BASIS, KERNEL, RESULT)                     \
    template BoundaryOperator<BASIS, RESULT>                            \
    modifiedHelmholtz3dSyntheticBoundaryOperator(                       \
        BoundaryOperator<BASIS, RESULT>                                 \
        (*)(                                                            \
            const shared_ptr<const Context<BASIS, RESULT> >& /*context*/, \
            const shared_ptr<const Space<BASIS> >& /*domain*/,          \
            const shared_ptr<const Space<BASIS> >& /*range*/,           \
            const shared_ptr<const Space<BASIS> >& /*dualToRange*/,     \
            KERNEL /*waveNumber*/,                                      \
            const std::string& /*label*/,                               \
            int /*symmetry*/,                                           \
            bool /*useInterpolation*/,                                  \
            int /*interpPtsPerWavelength*/),                            \
        const shared_ptr<const Context<BASIS, RESULT> >&,               \
        const shared_ptr<const Space<BASIS> >&,                         \
        const shared_ptr<const Space<BASIS> >&,                         \
        const shared_ptr<const Space<BASIS> >&,                         \
        KERNEL, std::string, int, bool, int, int)
FIBER_ITERATE_OVER_BASIS_KERNEL_AND_RESULT_TYPES(INSTANTIATE_FUNCTION);

} // namespace Bempp
