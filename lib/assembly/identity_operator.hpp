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

#include <boost/scoped_ptr.hpp>

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
    IdentityOperator(const IdentityOperator& other);
    virtual ~IdentityOperator() ;

    virtual int trialComponentCount() const { return 1; }

    virtual int testComponentCount() const { return 1; }

    virtual bool supportsRepresentation(AssemblyOptions::Representation repr) const;

private:
    virtual std::auto_ptr<LocalAssembler> makeAssemblerImpl(
            const LocalAssemblerFactory& assemblerFactory,
            const shared_ptr<const GeometryFactory>& testGeometryFactory,
            const shared_ptr<const GeometryFactory>& trialGeometryFactory,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
            const shared_ptr<const Fiber::OpenClHandler>& openClHandler,
            const ParallelisationOptions& parallelisationOptions,
            bool cacheSingularIntegrals) const;

    virtual std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleDetachedWeakFormImpl(
            const LocalAssemblerFactory& factory,
            const AssemblyOptions& options,
            Symmetry symmetry) const;

    virtual std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleDetachedWeakFormInternalImpl(
            LocalAssembler& assembler,
            const AssemblyOptions& options,
            Symmetry symmetry) const;

    std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleDetachedWeakFormInDenseMode(
            LocalAssembler& assembler,
            const AssemblyOptions& options,
            Symmetry symmetry) const;

    std::auto_ptr<DiscreteLinearOperator<ResultType> >
    assembleDetachedWeakFormInSparseMode(
            LocalAssembler& assembler,
            const AssemblyOptions& options,
            Symmetry symmetry) const;

private:
    struct Impl;
    boost::scoped_ptr<Impl> m_impl;
};

} // namespace Bempp

#endif
