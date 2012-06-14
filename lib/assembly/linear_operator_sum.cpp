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

#include "linear_operator_sum.hpp"

#include "discrete_linear_operator_sum.hpp"
#include "discrete_dense_linear_operator.hpp"

#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
LinearOperatorSum<BasisFunctionType, ResultType>::
LinearOperatorSum(const Base& term1, const Base& term2) :
    Base(term1.domain(), term1.range(), term1.dualToRange(),
         term1.label() + " + " + term2.label()),
    m_term1(term1.clone()), m_term2(term2.clone())
{
    assert(m_term1.get());
    assert(m_term2.get());
}

template <typename BasisFunctionType, typename ResultType>
LinearOperatorSum<BasisFunctionType, ResultType>::
LinearOperatorSum(const LinearOperatorSum& other) :
    Base(other), m_term1(other.m_term1->clone()), m_term2(other.m_term2->clone())
{
    assert(m_term1.get());
    assert(m_term2.get());
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<LinearOperator<BasisFunctionType, ResultType> >
LinearOperatorSum<BasisFunctionType, ResultType>::clone() const
{
    typedef LinearOperator<BasisFunctionType, ResultType> LinOp;
    typedef LinearOperatorSum<BasisFunctionType, ResultType> This;
    return std::auto_ptr<LinOp>(new This(*this));
}

template <typename BasisFunctionType, typename ResultType>
bool LinearOperatorSum<BasisFunctionType, ResultType>::supportsRepresentation(
        AssemblyOptions::Representation repr) const
{
    return (m_term1->supportsRepresentation(repr) &&
            m_term2->supportsRepresentation(repr));
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
LinearOperatorSum<BasisFunctionType, ResultType>::
assembleDetachedWeakFormImpl(const LocalAssemblerFactory& factory,
                             const AssemblyOptions& options,
                             Symmetry symmetry) const
{
    if (options.operatorRepresentation() == AssemblyOptions::DENSE)
        return assembleDetachedWeakFormInDenseMode(factory, options, symmetry);
    else
        return assembleDetachedWeakFormInArbitraryMode(factory, options, symmetry);
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
LinearOperatorSum<BasisFunctionType, ResultType>::
assembleDetachedWeakFormInDenseMode(
        const LocalAssemblerFactory& factory,
        const AssemblyOptions& options,
        Symmetry symmetry) const
{
    typedef DiscreteLinearOperator<ResultType> DiscreteLinOp;
    typedef DiscreteDenseLinearOperator<ResultType> DiscreteDenseLinOp;

    // The reset() calls below are made to release memory as early as possible
    std::auto_ptr<DiscreteLinOp> discreteTerm1 =
                m_term1->assembleDetachedWeakForm(factory, options, symmetry);
    arma::Mat<ResultType> sum = discreteTerm1->asMatrix();
    discreteTerm1.reset();
    std::auto_ptr<DiscreteLinOp> discreteTerm2 =
                m_term1->assembleDetachedWeakForm(factory, options, symmetry);
    sum += discreteTerm2->asMatrix();
    discreteTerm2.reset();

    return std::auto_ptr<DiscreteLinOp>(new DiscreteDenseLinOp(sum));
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<DiscreteLinearOperator<ResultType> >
LinearOperatorSum<BasisFunctionType, ResultType>::
assembleDetachedWeakFormInArbitraryMode(const LocalAssemblerFactory& factory,
                                        const AssemblyOptions& options,
                                        Symmetry symmetry) const
{
    typedef std::auto_ptr<DiscreteLinearOperator<ResultType> > DiscreteLinOpPtr;
    DiscreteLinOpPtr discreteTerm1 = 
        m_term1->assembleDetachedWeakForm(factory, options, symmetry);
    DiscreteLinOpPtr discreteTerm2 = 
        m_term2->assembleDetachedWeakForm(factory, options, symmetry);
    assert(discreteTerm1.get());
    assert(discreteTerm2.get());
    return DiscreteLinOpPtr(new DiscreteLinearOperatorSum<ResultType>(
                                discreteTerm1, discreteTerm2));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(LinearOperatorSum);

} // namespace Bempp
