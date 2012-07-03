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

#include "abstract_boundary_operator_sum.hpp"

#include "discrete_boundary_operator_sum.hpp"

#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
AbstractBoundaryOperatorSum<BasisFunctionType, ResultType>::
AbstractBoundaryOperatorSum(const Base& term1, const Base& term2) :
    Base(term1.domain(), term1.range(), term1.dualToRange(),
         term1.label() + " + " + term2.label()),
    m_term1(term1.clone()), m_term2(term2.clone())
{
    assert(m_term1.get());
    assert(m_term2.get());

    if (&m_term1->domain() != &m_term2->domain())
        throw std::invalid_argument(
                "AbstractBoundaryOperatorSum::AbstractBoundaryOperatorSum(" +
                m_term1->label() +
                ", " +
                m_term2->label() +
                "): Domains of the two terms must be equal");
    if (&m_term1->range() != &m_term2->range())
        throw std::invalid_argument(
                "AbstractBoundaryOperatorSum::AbstractBoundaryOperatorSum(" +
                m_term1->label() +
                ", " +
                m_term2->label() +
                "): Ranges of the two terms must be equal");
    if (&m_term1->dualToRange() != &m_term2->dualToRange())
        throw std::invalid_argument(
                "AbstractBoundaryOperatorSum::AbstractBoundaryOperatorSum(" +
                m_term1->label() +
                ", " +
                m_term2->label() +
                "): Spaces dual to the ranges of the two terms must be equal");
}

template <typename BasisFunctionType, typename ResultType>
AbstractBoundaryOperatorSum<BasisFunctionType, ResultType>::
AbstractBoundaryOperatorSum(const AbstractBoundaryOperatorSum& other) :
    Base(other), m_term1(other.m_term1->clone()), m_term2(other.m_term2->clone())
{
    assert(m_term1.get());
    assert(m_term2.get());
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<AbstractBoundaryOperator<BasisFunctionType, ResultType> >
AbstractBoundaryOperatorSum<BasisFunctionType, ResultType>::clone() const
{
    typedef AbstractBoundaryOperator<BasisFunctionType, ResultType> LinOp;
    typedef AbstractBoundaryOperatorSum<BasisFunctionType, ResultType> This;
    return std::auto_ptr<LinOp>(new This(*this));
}

template <typename BasisFunctionType, typename ResultType>
bool AbstractBoundaryOperatorSum<BasisFunctionType, ResultType>::supportsRepresentation(
        AssemblyOptions::Representation repr) const
{
    return (m_term1->supportsRepresentation(repr) &&
            m_term2->supportsRepresentation(repr));
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
AbstractBoundaryOperatorSum<BasisFunctionType, ResultType>::
assembleWeakFormImpl(const LocalAssemblerFactory& factory,
                     const AssemblyOptions& options,
                     Symmetry symmetry)
{
    typedef DiscreteBoundaryOperator<ResultType> DiscreteLinOp;

    if (!m_term1->isWeakFormAssembled())
        m_term1->assembleWeakForm(factory, options, symmetry);
    shared_ptr<const DiscreteLinOp> discreteTerm1 = m_term1->weakForm();
    assert(discreteTerm1);

    if (!m_term2->isWeakFormAssembled())
        m_term2->assembleWeakForm(factory, options, symmetry);
    shared_ptr<const DiscreteLinOp> discreteTerm2 = m_term2->weakForm();
    assert(discreteTerm2);

    return shared_ptr<DiscreteLinOp>(new DiscreteBoundaryOperatorSum<ResultType>(
                                         discreteTerm1, discreteTerm2));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AbstractBoundaryOperatorSum);

} // namespace Bempp
