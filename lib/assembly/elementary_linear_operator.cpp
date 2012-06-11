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

#include "elementary_linear_operator.hpp"

#include "discrete_linear_operator.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/local_assembler_for_operators.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
ElementaryLinearOperator<BasisFunctionType, ResultType>::
ElementaryLinearOperator(const Space<BasisFunctionType>& testSpace,
                         const Space<BasisFunctionType>& trialSpace) :
    LinearOperator<BasisFunctionType, ResultType>(testSpace, trialSpace)
{
}

template <typename BasisFunctionType, typename ResultType>
ElementaryLinearOperator<BasisFunctionType, ResultType>::
~ElementaryLinearOperator()
{
}

template <typename BasisFunctionType, typename ResultType>
std::vector<const ElementaryLinearOperator<BasisFunctionType, ResultType>*>
ElementaryLinearOperator<BasisFunctionType, ResultType>::
constituentOperators() const
{
    return std::vector<const ElementaryLinearOperator<
            BasisFunctionType, ResultType>*>(1 /* size */, this);
}

template <typename BasisFunctionType, typename ResultType>
std::vector<ResultType>
ElementaryLinearOperator<BasisFunctionType, ResultType>::
constituentOperatorWeights() const
{
    return std::vector<ResultType>(1 /* size */, 1.);
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
ElementaryLinearOperator<BasisFunctionType, ResultType>::
assembleDetachedWeakFormInternal(
        LocalAssembler& assembler,
        const AssemblyOptions& options,
        Symmetry symmetry) const
{
    return assembleDetachedWeakFormInternalImpl(assembler, options, symmetry);
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<typename ElementaryLinearOperator<BasisFunctionType, ResultType>::LocalAssembler>
ElementaryLinearOperator<BasisFunctionType, ResultType>::makeAssembler(
        const LocalAssemblerFactory& assemblerFactory,
        const shared_ptr<const GeometryFactory>& testGeometryFactory,
        const shared_ptr<const GeometryFactory>& trialGeometryFactory,
        const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
        const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
        const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
        const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
        const shared_ptr<const Fiber::OpenClHandler>& openClHandler,
        const ParallelisationOptions& parallelisationOptions,
        bool cacheSingularIntegrals) const
{
    return makeAssemblerImpl(assemblerFactory,
                             testGeometryFactory, trialGeometryFactory,
                             testRawGeometry, trialRawGeometry,
                             testBases, trialBases, openClHandler,
                             parallelisationOptions,
                             cacheSingularIntegrals);
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<typename ElementaryLinearOperator<BasisFunctionType, ResultType>::LocalAssembler>
ElementaryLinearOperator<BasisFunctionType, ResultType>::makeAssembler(
        const LocalAssemblerFactory& assemblerFactory,
        const AssemblyOptions& options) const
{
    typedef Fiber::RawGridGeometry<CoordinateType> RawGridGeometry;
    typedef std::vector<const Fiber::Basis<BasisFunctionType>*> BasisPtrVector;

    shared_ptr<RawGridGeometry> testRawGeometry, trialRawGeometry;
    shared_ptr<GeometryFactory> testGeometryFactory, trialGeometryFactory;
    shared_ptr<Fiber::OpenClHandler> openClHandler;
    shared_ptr<BasisPtrVector> testBases, trialBases;
    bool cacheSingularIntegrals;

    std::cout << "Collecting data for assembler construction" << std::endl;
    collectDataForAssemblerConstruction(options,
                                        testRawGeometry, trialRawGeometry,
                                        testGeometryFactory, trialGeometryFactory,
                                        testBases, trialBases,
                                        openClHandler, cacheSingularIntegrals);
    std::cout << "Collection finished." << std::endl;

    return makeAssemblerImpl(assemblerFactory,
                             testGeometryFactory, trialGeometryFactory,
                             testRawGeometry, trialRawGeometry,
                             testBases, trialBases, openClHandler,
                             options.parallelisationOptions(),
                             cacheSingularIntegrals);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(ElementaryLinearOperator);

} // namespace Bempp
