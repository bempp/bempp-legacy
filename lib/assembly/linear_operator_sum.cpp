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
shared_ptr<DiscreteLinearOperator<ResultType> >
LinearOperatorSum<BasisFunctionType, ResultType>::
assembleWeakFormImpl(const LocalAssemblerFactory& factory,
                     const AssemblyOptions& options,
                     Symmetry symmetry)
{
    typedef DiscreteLinearOperator<ResultType> DiscreteLinOp;

    if (!m_term1->isWeakFormAssembled())
        m_term1->assembleWeakForm(factory, options, symmetry);
    shared_ptr<const DiscreteLinOp> discreteTerm1 = m_term1->weakForm();
    assert(discreteTerm1);

    if (!m_term2->isWeakFormAssembled())
        m_term2->assembleWeakForm(factory, options, symmetry);
    shared_ptr<const DiscreteLinOp> discreteTerm2 = m_term2->weakForm();
    assert(discreteTerm2);

    return shared_ptr<DiscreteLinOp>(new DiscreteLinearOperatorSum<ResultType>(
                                         discreteTerm1, discreteTerm2));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(LinearOperatorSum);

} // namespace Bempp
