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

#include "../common/config_trilinos.hpp"

#include "boundary_operator_composition.hpp"

#include "discrete_boundary_operator_composition.hpp"
#include "mass_matrix_container.hpp"
#include "mass_matrix_container_initialiser.hpp"

#include "../common/shared_ptr.hpp"
#include "../fiber/explicit_instantiation.hpp"

#ifdef WITH_TRILINOS

#include <boost/make_shared.hpp>
#endif

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
BoundaryOperatorComposition<BasisFunctionType, ResultType>::
BoundaryOperatorComposition(const Base& outer, const Base& inner) :
    Base(outer.domain(), outer.range(), outer.dualToRange(),
         outer.label() + " * " + inner.label()),
    m_outer(outer.clone()), m_inner(inner.clone())
{
    assert(m_outer.get());
    assert(m_inner.get());

    if (&m_outer->domain() != &m_inner->range())
        throw std::invalid_argument(
                "BoundaryOperatorComposition::BoundaryOperatorComposition(" +
                m_outer->label() +
                ", " +
                m_inner->label() +
                "): Domain of the outer term must be equal to the range of "
                "the inner term");
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperatorComposition<BasisFunctionType, ResultType>::
BoundaryOperatorComposition(const BoundaryOperatorComposition& other) :
    Base(other), m_outer(other.m_outer->clone()), m_inner(other.m_inner->clone())
{
    assert(m_outer.get());
    assert(m_inner.get());
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<BoundaryOperator<BasisFunctionType, ResultType> >
BoundaryOperatorComposition<BasisFunctionType, ResultType>::clone() const
{
    typedef BoundaryOperator<BasisFunctionType, ResultType> LinOp;
    typedef BoundaryOperatorComposition<BasisFunctionType, ResultType> This;
    return std::auto_ptr<LinOp>(new This(*this));
}

template <typename BasisFunctionType, typename ResultType>
bool BoundaryOperatorComposition<BasisFunctionType, ResultType>::supportsRepresentation(
        AssemblyOptions::Representation repr) const
{
    return (m_outer->supportsRepresentation(repr) &&
            m_inner->supportsRepresentation(repr));
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
BoundaryOperatorComposition<BasisFunctionType, ResultType>::
assembleWeakFormImpl(const LocalAssemblerFactory& factory,
                     const AssemblyOptions& options,
                     Symmetry symmetry)
{
    typedef DiscreteBoundaryOperator<ResultType> DiscreteLinOp;

    if (!m_outer->isWeakFormAssembled())
        m_outer->assembleWeakForm(factory, options, symmetry);
    shared_ptr<const DiscreteLinOp> discreteOuter = m_outer->weakForm();
    assert(discreteOuter);

    if (!m_inner->isWeakFormAssembled())
        m_inner->assembleWeakForm(factory, options, symmetry);
    shared_ptr<const DiscreteLinOp> discreteInner = m_inner->weakForm();
    assert(discreteInner);

    // Calculate the (pseudo)inverse mass matrix
    MassMatrixContainerInitialiser<BasisFunctionType, ResultType> mmcInitialiser(
                m_inner->range(), m_inner->dualToRange());
    std::auto_ptr<MassMatrixContainer<ResultType> > mmc = mmcInitialiser();

    shared_ptr<const DiscreteLinOp> temp =
            boost::make_shared<DiscreteBoundaryOperatorComposition<ResultType> >(
                mmc->massMatrixPseudoinverse, discreteInner);
    return boost::make_shared<DiscreteBoundaryOperatorComposition<ResultType> >(
                discreteOuter, temp);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(BoundaryOperatorComposition);

} // namespace Bempp
