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

#include "null_operator.hpp"

#include "boundary_operator.hpp"
#include "discrete_null_boundary_operator.hpp"
#include "context.hpp"

#include "../common/boost_make_shared_fwd.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <stdexcept>

namespace Bempp
{

////////////////////////////////////////////////////////////////////////////////
// NullOperatorId

template <typename BasisFunctionType, typename ResultType>
NullOperatorId<BasisFunctionType, ResultType>::NullOperatorId(
        const NullOperator<BasisFunctionType, ResultType>& op) :
    m_domain(op.domain().get()), m_range(op.range().get()),
    m_dualToRange(op.dualToRange().get())
{
}

template <typename BasisFunctionType, typename ResultType>
size_t NullOperatorId<BasisFunctionType, ResultType>::hash() const
{
    typedef NullOperator<BasisFunctionType, ResultType>
            OperatorType;
    size_t result = tbb::tbb_hasher(typeid(OperatorType).name());
    tbb_hash_combine(result, m_domain);
    tbb_hash_combine(result, m_range);
    tbb_hash_combine(result, m_dualToRange);
    return result;
}

template <typename BasisFunctionType, typename ResultType>
bool NullOperatorId<BasisFunctionType, ResultType>::isEqual(
        const AbstractBoundaryOperatorId &other) const
{
    // dynamic_cast won't suffice since we want to make sure both objects
    // are of exactly the same type (dynamic_cast would succeed for a subclass)
    if (typeid(other) == typeid(*this))
    {
        const NullOperatorId& otherCompatible =
            static_cast<const NullOperatorId&>(other);
        return (m_domain == otherCompatible.m_domain &&
                m_range == otherCompatible.m_range &&
                m_dualToRange == otherCompatible.m_dualToRange);
    }
    else
        return false;
}

template <typename BasisFunctionType, typename ResultType>
NullOperator<BasisFunctionType, ResultType>::NullOperator(
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label,
        int symmetry) :
    Base(domain, range, dualToRange, label,
         (symmetry & AUTO_SYMMETRY ?
              (domain == dualToRange ?
                   SYMMETRIC | HERMITIAN :
                   NO_SYMMETRY) :
              symmetry)),
    m_id(boost::make_shared<NullOperatorId<BasisFunctionType, ResultType> >(
             *this))
{
}

template <typename BasisFunctionType, typename ResultType>
NullOperator<BasisFunctionType, ResultType>::~NullOperator()
{
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const AbstractBoundaryOperatorId>
NullOperator<BasisFunctionType, ResultType>::id() const
{
    return m_id;
}

template <typename BasisFunctionType, typename ResultType>
bool NullOperator<BasisFunctionType, ResultType>::isLocal() const
{
    return true;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
NullOperator<BasisFunctionType, ResultType>::assembleWeakFormImpl(
        const Context<BasisFunctionType, ResultType>& context) const
{
    return reallyAssembleWeakForm();
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<DiscreteBoundaryOperator<ResultType> >
NullOperator<BasisFunctionType, ResultType>::reallyAssembleWeakForm() const
{
    shared_ptr<DiscreteBoundaryOperator<ResultType> > result(
                new DiscreteNullBoundaryOperator<ResultType>(
                    this->dualToRange()->globalDofCount(),
                    this->domain()->globalDofCount()));
    return result;
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
nullOperator(const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
             const shared_ptr<const Space<BasisFunctionType> >& domain,
             const shared_ptr<const Space<BasisFunctionType> >& range,
             const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
             const std::string& label,
             int symmetry)
{
    typedef NullOperator<BasisFunctionType, ResultType> NullOp;
    return BoundaryOperator<BasisFunctionType, ResultType>(
                context,
                boost::make_shared<NullOp>(domain, range, dualToRange,
                                       label, symmetry));
}

#define INSTANTIATE_NONMEMBER_CONSTRUCTOR(BASIS, RESULT) \
    template BoundaryOperator<BASIS, RESULT> \
    nullOperator( \
        const shared_ptr<const Context<BASIS, RESULT> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const shared_ptr<const Space<BASIS> >&, \
        const std::string&, \
        int)
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_NONMEMBER_CONSTRUCTOR);

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(NullOperator);

} // namespace Bempp
