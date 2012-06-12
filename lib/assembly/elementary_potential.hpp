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

#ifndef bempp_elementary_potential_hpp
#define bempp_elementary_potential_hpp

#include "potential.hpp"

namespace Fiber
{

template <typename CoordinateType> class CollectionOfBasisTransformations;
template <typename KernelType> class CollectionOfKernels;
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class KernelTrialIntegral;
template <typename ResultType> class EvaluatorForIntegralOperators;

} // namespace Bempp

namespace Bempp
{

template <typename BasisFunctionType_, typename KernelType_, typename ResultType_>
class ElementaryPotential : public Potential<BasisFunctionType_, ResultType_>
{
    typedef Potential<BasisFunctionType_, ResultType_> Base;
public:
    typedef typename Base::BasisFunctionType BasisFunctionType;
    typedef typename Base::ResultType ResultType;
    typedef KernelType_ KernelType;
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::LocalAssemblerFactory LocalAssemblerFactory;
    typedef Fiber::EvaluatorForIntegralOperators<ResultType> Evaluator;
    typedef Fiber::CollectionOfBasisTransformations<CoordinateType>
    CollectionOfBasisTransformations;
    typedef Fiber::CollectionOfKernels<KernelType> CollectionOfKernels;
    typedef Fiber::KernelTrialIntegral<BasisFunctionType, KernelType, ResultType>
    KernelTrialIntegral;

    virtual std::auto_ptr<InterpolatedFunction<ResultType_> > evaluateOnGrid(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const Grid& evaluationGrid,
            const LocalAssemblerFactory& assemblerFactory,
            const EvaluationOptions& options) const;

    virtual arma::Mat<ResultType_> evaluateAtPoints(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const arma::Mat<CoordinateType>& evaluationPoints,
            const LocalAssemblerFactory& assemblerFactory,
            const EvaluationOptions& options) const;

private:
    std::auto_ptr<Evaluator>
    makeEvaluator(
            const GridFunction<BasisFunctionType, ResultType>& argument,
            const LocalAssemblerFactory& factory,
            const EvaluationOptions& options) const;

    virtual const CollectionOfKernels& kernels() const = 0;
    virtual const CollectionOfBasisTransformations&
    trialTransformations() const = 0;
    virtual const KernelTrialIntegral& integral() const = 0;
};

} // namespace Bempp

#endif
