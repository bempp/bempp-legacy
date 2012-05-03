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


#ifndef bempp_identity_operator_hpp
#define bempp_identity_operator_hpp

#include "elementary_linear_operator.hpp"
#include "../fiber/scalar_function_value.hpp"

namespace Fiber
{

template <typename ResultType> class LocalAssemblerForOperators;

} // namespace Fiber

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
class IdentityOperator : public ElementaryLinearOperator<BasisFunctionType, ResultType>
{
    typedef ElementaryLinearOperator<BasisFunctionType, ResultType> Base;
public:
    typedef typename Base::LocalAssemblerFactory LocalAssemblerFactory;
    typedef typename Base::LocalAssembler LocalAssembler;
    typedef typename Base::CoordinateType CoordinateType;

    IdentityOperator(const Space<BasisFunctionType>& testSpace,
                     const Space<BasisFunctionType>& trialSpace);

    virtual int trialComponentCount() const { return 1; }

    virtual int testComponentCount() const { return 1; }

    virtual bool supportsRepresentation(AssemblyOptions::Representation repr) const;

    virtual std::auto_ptr<LocalAssembler> makeAssembler(
            const LocalAssemblerFactory& assemblerFactory,
            const GeometryFactory& geometryFactory,
            const Fiber::RawGridGeometry<CoordinateType>& rawGeometry,
            const std::vector<const Fiber::Basis<BasisFunctionType>*>& testBases,
            const std::vector<const Fiber::Basis<BasisFunctionType>*>& trialBases,
            const Fiber::OpenClHandler& openClHandler,
            bool cacheSingularIntegrals) const;

    virtual std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakForm(
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options) const;

    virtual std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakFormInternal(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

private:
    std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakFormInDenseMode(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

    std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleWeakFormInSparseMode(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

private:
    Fiber::ScalarFunctionValue<CoordinateType> m_expression;
};

} // namespace Bempp

#endif
