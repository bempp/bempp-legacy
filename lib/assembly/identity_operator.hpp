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

#include "../common/common.hpp"

#include "elementary_abstract_boundary_operator.hpp"

#include "abstract_boundary_operator_id.hpp"
#include <boost/scoped_ptr.hpp>

namespace Fiber
{

template <typename ResultType> class LocalAssemblerForOperators;

} // namespace Fiber

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType> class BoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class IdentityOperator;

template <typename BasisFunctionType, typename ResultType>
class IdentityOperatorId : public AbstractBoundaryOperatorId
{
public:
    IdentityOperatorId(const IdentityOperator<BasisFunctionType, ResultType>& op);
    virtual size_t hash() const;
    virtual bool isEqual(const AbstractBoundaryOperatorId &other) const;

private:
    const Space<BasisFunctionType>* m_domain;
    const Space<BasisFunctionType>* m_range;
    const Space<BasisFunctionType>* m_dualToRange;
};

template <typename BasisFunctionType_, typename ResultType_>
class IdentityOperator :
        public ElementaryAbstractBoundaryOperator<BasisFunctionType_, ResultType_>
{
    typedef ElementaryAbstractBoundaryOperator<BasisFunctionType_, ResultType_> Base;
public:
    typedef typename Base::BasisFunctionType BasisFunctionType;
    typedef typename Base::ResultType ResultType;
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::QuadratureStrategy QuadratureStrategy;
    typedef typename Base::LocalAssembler LocalAssembler;

    IdentityOperator(const shared_ptr<const Space<BasisFunctionType> >& domain,
                     const shared_ptr<const Space<BasisFunctionType> >& range,
                     const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
                     const std::string& label = "");
    IdentityOperator(const IdentityOperator& other);
    virtual ~IdentityOperator();

    virtual shared_ptr<const AbstractBoundaryOperatorId> id() const;

    virtual bool supportsRepresentation(AssemblyOptions::Representation repr) const;

protected:
    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormImpl(
            const Context<BasisFunctionType, ResultType>& context) const;

private:
    virtual std::auto_ptr<LocalAssembler> makeAssemblerImpl(
            const QuadratureStrategy& quadStrategy,
            const shared_ptr<const GeometryFactory>& testGeometryFactory,
            const shared_ptr<const GeometryFactory>& trialGeometryFactory,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& testRawGeometry,
            const shared_ptr<const Fiber::RawGridGeometry<CoordinateType> >& trialRawGeometry,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& testBases,
            const shared_ptr<const std::vector<const Fiber::Basis<BasisFunctionType>*> >& trialBases,
            const shared_ptr<const Fiber::OpenClHandler>& openClHandler,
            const ParallelizationOptions& parallelizationOptions,
            bool cacheSingularIntegrals) const;

    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInternalImpl(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

    std::auto_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInDenseMode(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

    std::auto_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormInSparseMode(
            LocalAssembler& assembler,
            const AssemblyOptions& options) const;

private:
    struct Impl;
    boost::scoped_ptr<Impl> m_impl;
    shared_ptr<const AbstractBoundaryOperatorId> m_id;
};

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
identityOperator(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
                 const shared_ptr<const Space<BasisFunctionType> >& domain,
                 const shared_ptr<const Space<BasisFunctionType> >& range,
                 const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
                 const std::string& label = "");

} // namespace Bempp

#endif
