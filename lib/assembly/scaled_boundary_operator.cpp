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

#include "scaled_boundary_operator.hpp"

#include "scaled_discrete_boundary_operator.hpp"

#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
ScaledBoundaryOperator<BasisFunctionType, ResultType>::
ScaledBoundaryOperator(ResultType weight, const Base& linOp) :
    Base(linOp.domain(), linOp.range(), linOp.dualToRange(),
         toString(weight) + " * " + linOp.label()),
    m_weight(weight), m_operator(linOp.clone())
{
    assert(m_operator.get());
}

template <typename BasisFunctionType, typename ResultType>
ScaledBoundaryOperator<BasisFunctionType, ResultType>::
ScaledBoundaryOperator(const ScaledBoundaryOperator& other) :
    Base(other), m_weight(other.m_weight), m_operator(other.m_operator->clone())
{
    assert(m_operator.get());
}

template <typename BasisFunctionType, typename ResultType>
std::auto_ptr<BoundaryOperator<BasisFunctionType, ResultType> >
ScaledBoundaryOperator<BasisFunctionType, ResultType>::clone() const
{
    typedef BoundaryOperator<BasisFunctionType, ResultType> LinOp;
    typedef ScaledBoundaryOperator<BasisFunctionType, ResultType> This;
    return std::auto_ptr<LinOp>(new This(*this));
}

template <typename BasisFunctionType, typename ResultType>
bool ScaledBoundaryOperator<BasisFunctionType, ResultType>::supportsRepresentation(
        AssemblyOptions::Representation repr) const
{
    return m_operator->supportsRepresentation(repr);
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
ScaledBoundaryOperator<BasisFunctionType, ResultType>::
assembleWeakFormImpl(const LocalAssemblerFactory& factory,
                     const AssemblyOptions& options,
                     Symmetry symmetry)
{
    if (!m_operator->isWeakFormAssembled())
        m_operator->assembleWeakForm(factory, options, symmetry);
    shared_ptr<const DiscreteBoundaryOperator<ResultType> > wrappedDiscreteLinOp =
            m_operator->weakForm();
    return shared_ptr<DiscreteBoundaryOperator<ResultType> >(
                new ScaledDiscreteBoundaryOperator<ResultType>(
                    m_weight, wrappedDiscreteLinOp));
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(ScaledBoundaryOperator);

} // namespace Bempp
