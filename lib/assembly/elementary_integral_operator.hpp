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

#ifndef bempp_elementary_integral_operator_hpp
#define bempp_elementary_integral_operator_hpp

#include "elementary_linear_operator.hpp"
#include "../common/multidimensional_arrays.hpp"
#include "../common/types.hpp"
#include "../fiber/kernel.hpp"
#include "../fiber/types.hpp"

#include <vector>
#include <armadillo>

namespace Fiber
{

template <typename CoordinateType> class Expression;
template <typename ResultType> class LocalAssemblerForOperators;
template <typename ResultType> class EvaluatorForIntegralOperators;

} // namespace Fiber

namespace Bempp
{

class EvaluationOptions;
template <typename BasisFunctionType, typename ResultType> class GridFunction;
template <typename ValueType> class InterpolatedFunction;
template <typename BasisFunctionType, typename ResultType> class WeakFormAcaAssemblyHelper;

template <typename BasisFunctionType, typename KernelType, typename ResultType>
class ElementaryIntegralOperator :
        public ElementaryLinearOperator<BasisFunctionType, ResultType>
{
    typedef ElementaryLinearOperator<BasisFunctionType, ResultType> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::LocalAssemblerFactory LocalAssemblerFactory;
    typedef typename Base::LocalAssembler LocalAssembler;
    typedef Fiber::EvaluatorForIntegralOperators<ResultType> Evaluator;

    ElementaryIntegralOperator(const Space<BasisFunctionType> &testSpace,
                               const Space<BasisFunctionType> &trialSpace);

    virtual int trialComponentCount() const {
        return kernel().domainDimension();
    }

    virtual int testComponentCount() const {
        return kernel().codomainDimension();
    }

    virtual bool supportsRepresentation(AssemblyOptions::Representation repr) const;

    virtual bool isRegular() const = 0;

    virtual std::auto_ptr<LocalAssembler> makeAssembler(
            const LocalAssemblerFactory& assemblerFactory,
            const GeometryFactory& geometryFactory,
            const Fiber::RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Fiber::Basis<BasisFunctionType>*>& testBases,
            const std::vector<const Fiber::Basis<BasisFunctionType>*>& trialBases,
            const Fiber::OpenClHandler& openClHandler,
            const ParallelisationOptions& parallelisationOptions,
            bool cacheSingularIntegrals) const;

    virtual std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakForm(
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options) const;

    virtual std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakFormInternal(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

    // We might define a superclass IntegralOperator that might represent
    // a superposition of elementary linear operators (defined at points
    // off surface). Then the virtual attribute here would be useful.
    virtual std::auto_ptr<InterpolatedFunction<ResultType> > applyOffSurface(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const Grid& evaluationGrid,
            const LocalAssemblerFactory& factory,
            const EvaluationOptions& options) const;

    // TODO: define applyOnSurface() for *all* operators (including Id).

private:
    virtual const Fiber::Kernel<KernelType>& kernel() const = 0;
    virtual const Fiber::Expression<CoordinateType>& testExpression() const = 0;
    virtual const Fiber::Expression<CoordinateType>& trialExpression() const = 0;

    /** @}
        \name Weak form assembly
        @{ */
    std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakFormInDenseMode(
            LocalAssembler& assembler,
            const AssemblyOptions &options) const;
    std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakFormInAcaMode(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;
    /** @} */

    std::auto_ptr<InterpolatedFunction<ResultType> >
    applyOffSurfaceWithKnownEvaluator(
            const Grid& evaluationGrid,
            const Evaluator& evaluator,
            const EvaluationOptions& options) const;
};

} // namespace Bempp

#endif
