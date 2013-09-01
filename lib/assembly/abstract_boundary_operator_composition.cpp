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

#include "bempp/common/config_trilinos.hpp"

#include "abstract_boundary_operator_composition.hpp"

#include "abstract_boundary_operator_pseudoinverse.hpp"
#include "boundary_operator.hpp"
#include "context.hpp"
#include "discrete_boundary_operator_composition.hpp"
#include "identity_operator.hpp"

#include "../common/boost_make_shared_fwd.hpp"
#include "../common/shared_ptr.hpp"
#include "../fiber/explicit_instantiation.hpp"

#ifdef WITH_TRILINOS

#include <boost/make_shared.hpp>
#endif

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
AbstractBoundaryOperatorComposition<BasisFunctionType, ResultType>::
AbstractBoundaryOperatorComposition(
        const BoundaryOperator<BasisFunctionType, ResultType>& outer,
        const BoundaryOperator<BasisFunctionType, ResultType>& inner,
        int symmetry) :
    Base(inner.domain(), outer.range(), outer.dualToRange(),
         "(" + outer.label() + ") * (" + inner.label() + ")",
         symmetry),
    m_outer(outer), m_inner(inner)
{
    assert(m_outer.abstractOperator());
    assert(m_inner.abstractOperator());

    if (!(m_outer.domain()->spaceIsCompatible(*m_inner.range())))
        throw std::invalid_argument(
                "AbstractBoundaryOperatorComposition::AbstractBoundaryOperatorComposition(" +
                m_outer.label() +
                ", " +
                m_inner.label() +
                "): Domain of the outer term must be compatible to the range of "
                "the inner term");
}

template <typename BasisFunctionType, typename ResultType>
bool AbstractBoundaryOperatorComposition<BasisFunctionType, ResultType>::isLocal() const
{
    return (m_outer.abstractOperator()->isLocal() &&
            m_inner.abstractOperator()->isLocal());
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
AbstractBoundaryOperatorComposition<BasisFunctionType, ResultType>::
assembleWeakFormImpl(const Context<BasisFunctionType, ResultType>& context) const
{
    typedef BoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;
    typedef DiscreteBoundaryOperator<ResultType> DiscreteLinOp;

    shared_ptr<const DiscreteLinOp> discreteOuter = m_outer.weakForm();
    shared_ptr<const DiscreteLinOp> discreteInner = m_inner.weakForm();

    // Calculate the (pseudo)inverse mass matrix
    BoundaryOp id = identityOperator(
                // We don't need a persistent shared_ptr since identityOperator
                // will go out of scope at the end of this function anyway.
                // All we need is a weak form.
                make_shared_from_ref(context),
                m_inner.range(), m_inner.range(), m_inner.dualToRange());
    BoundaryOp pinvId = pseudoinverse(id,m_inner.dualToRange());
    // Dual space not important here. Could be anything.

    shared_ptr<const DiscreteLinOp> temp =
            boost::make_shared<DiscreteBoundaryOperatorComposition<ResultType> >(
                pinvId.weakForm(), discreteInner);
    return boost::make_shared<DiscreteBoundaryOperatorComposition<ResultType> >(
                discreteOuter, temp);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(AbstractBoundaryOperatorComposition);

} // namespace Bempp
