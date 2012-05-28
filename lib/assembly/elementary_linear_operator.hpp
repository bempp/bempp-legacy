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


#ifndef bempp_elementary_linear_operator_hpp
#define bempp_elementary_linear_operator_hpp

#include "linear_operator.hpp"

#include "../common/shared_ptr.hpp"
#include <stdexcept>
#include <vector>

namespace Fiber
{

template <typename ResultType> class LocalAssemblerForOperators;
template <typename CoordinateType> class RawGridGeometry;
template <typename ValueType> class Basis;
class OpenClHandler;

} // namespace Fiber

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
class ElementaryLinearOperator : public LinearOperator<BasisFunctionType, ResultType>
{
    typedef LinearOperator<BasisFunctionType, ResultType> Base;
public:
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::LocalAssemblerFactory LocalAssemblerFactory;
    typedef Fiber::LocalAssemblerForOperators<ResultType> LocalAssembler;

    ElementaryLinearOperator(const Space<BasisFunctionType>& testSpace,
                             const Space<BasisFunctionType>& trialSpace);

    /** \brief Using a specified factory, construct a local assembler suitable
       for this operator. */
    virtual std::auto_ptr<LocalAssembler> makeAssembler(
            const LocalAssemblerFactory& assemblerFactory,
            const shared_ptr<const GeometryFactory>& testGeometryFactory,
            const shared_ptr<const GeometryFactory>& trialGeometryFactory,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
            const shared_ptr<const Fiber::OpenClHandler>& openClHandler,
            const ParallelisationOptions& parallelisationOptions,
            bool cacheSingularIntegrals) const = 0;

    std::auto_ptr<LocalAssembler> makeAssemblerFromScratch(
            const LocalAssemblerFactory& assemblerFactory,
            const AssemblyOptions& options) const;

    /** \brief Assemble the operator's weak form using a specified local assembler.

      This function is not intended to be called directly by the user. */
    std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleDetachedWeakFormInternal(
            LocalAssembler& assembler,
            const AssemblyOptions& options,
            Symmetry symmetry = UNSYMMETRIC) const;

private:
    virtual std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleDetachedWeakFormInternalImpl(
            LocalAssembler& assembler,
            const AssemblyOptions& options,
            Symmetry symmetry) const = 0;
};

} // namespace Bempp

#endif
