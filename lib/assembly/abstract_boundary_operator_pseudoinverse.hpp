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

#ifndef bempp_abstract_operator_pseudoinverse_hpp
#define bempp_abstract_operator_pseudoinverse_hpp

#include "abstract_boundary_operator.hpp"
#include "abstract_boundary_operator_id.hpp"
#include "boundary_operator.hpp"

namespace Bempp
{

template <typename ResultType> class DiscreteDenseBoundaryOperator;
template <typename ResultType> class DiscreteSparseBoundaryOperator;

template <typename BasisFunctionType, typename ResultType>
class AbstractBoundaryOperatorPseudoinverseId : public AbstractBoundaryOperatorId
{
public:
    explicit AbstractBoundaryOperatorPseudoinverseId(
            const BoundaryOperator<BasisFunctionType, ResultType>& op);
    virtual size_t hash() const;
    virtual void dump() const;
    virtual bool isEqual(const AbstractBoundaryOperatorId &other) const;

private:
    shared_ptr<const AbstractBoundaryOperatorId> m_operatorToInvertId;
};

template <typename BasisFunctionType_, typename ResultType_>
class AbstractBoundaryOperatorPseudoinverse :
        public AbstractBoundaryOperator<BasisFunctionType_, ResultType_>
{
    typedef AbstractBoundaryOperator<BasisFunctionType_, ResultType_> Base;
public:
    /** \copydoc AbstractBoundaryOperator::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc AbstractBoundaryOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc AbstractBoundaryOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc AbstractBoundaryOperator::QuadratureStrategy */
    typedef typename Base::QuadratureStrategy QuadratureStrategy;

    explicit AbstractBoundaryOperatorPseudoinverse(
            const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOp);

    virtual bool isLocal() const;

    virtual shared_ptr<const AbstractBoundaryOperatorId> id() const;

    // TODO: implement id() (will consist of typeid(*this) and m_operator.abo()->id().

protected:
    virtual shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormImpl(
            const Context<BasisFunctionType, ResultType>& context) const;

private:
    shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormForSparseOperator(
            const Context<BasisFunctionType, ResultType>& context,
            const shared_ptr<const DiscreteSparseBoundaryOperator<ResultType> >&
            wrappedDiscreteOp) const;
    shared_ptr<DiscreteBoundaryOperator<ResultType_> >
    assembleWeakFormForDenseOperator(
            const Context<BasisFunctionType, ResultType>& context,
            const shared_ptr<const DiscreteDenseBoundaryOperator<ResultType> >&
            wrappedDiscreteOp) const;

private:
    BoundaryOperator<BasisFunctionType, ResultType> m_operator;
    shared_ptr<const AbstractBoundaryOperatorId> m_id;
};

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
pseudoinverse(const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOp);

} // namespace Bempp

#endif
