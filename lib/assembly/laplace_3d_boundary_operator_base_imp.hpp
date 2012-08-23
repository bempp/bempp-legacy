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

#ifndef bempp_laplace_3d_boundary_operator_base_imp_hpp
#define bempp_laplace_3d_boundary_operator_base_imp_hpp

#include "laplace_3d_boundary_operator_base.hpp"
#include "abstract_boundary_operator_id.hpp"
#include "../common/boost_make_shared_fwd.hpp"

namespace Bempp
{

////////////////////////////////////////////////////////////////////////////////
// Laplace3dBoundaryOperatorId

template <typename BasisFunctionType>
class Laplace3dBoundaryOperatorId : public AbstractBoundaryOperatorId
{
public:
    template <typename Impl, typename ResultType>
    explicit Laplace3dBoundaryOperatorId(
            const Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>& op);
    virtual size_t hash() const;
    virtual void dump() const;
    virtual bool isEqual(const AbstractBoundaryOperatorId &other) const;

private:
    const std::type_info& m_typeInfo;
    const Space<BasisFunctionType>* m_domain;
    const Space<BasisFunctionType>* m_range;
    const Space<BasisFunctionType>* m_dualToRange;
};

template <typename BasisFunctionType>
template <typename Impl, typename ResultType>
Laplace3dBoundaryOperatorId<BasisFunctionType>::Laplace3dBoundaryOperatorId(
        const Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>& op) :
    m_typeInfo(typeid(op)),
    m_domain(op.domain().get()), m_range(op.range().get()),
    m_dualToRange(op.dualToRange().get())
{
}

template <typename BasisFunctionType>
size_t Laplace3dBoundaryOperatorId<BasisFunctionType>::hash() const
{
    size_t result = tbb::tbb_hasher(m_typeInfo.name());
    tbb_hash_combine(result, m_domain);
    tbb_hash_combine(result, m_range);
    tbb_hash_combine(result, m_dualToRange);
    return result;
}

template <typename BasisFunctionType>
void Laplace3dBoundaryOperatorId<BasisFunctionType>::dump() const
{
    std::cout << m_typeInfo.name() << ", " << m_domain << ", "
              << m_range << ", " << m_dualToRange << std::endl;
}

template <typename BasisFunctionType>
bool Laplace3dBoundaryOperatorId<BasisFunctionType>::isEqual(
        const AbstractBoundaryOperatorId &other) const
{
    // dynamic_cast won't suffice since we want to make sure both objects
    // are of exactly the same type (dynamic_cast would succeed for a subclass)
    if (typeid(other) == typeid(*this)) {
        const Laplace3dBoundaryOperatorId& otherCompatible =
            static_cast<const Laplace3dBoundaryOperatorId&>(other);
        return (m_typeInfo == otherCompatible.m_typeInfo &&
                m_domain == otherCompatible.m_domain &&
                m_range == otherCompatible.m_range &&
                m_dualToRange == otherCompatible.m_dualToRange);
    }
    else
        return false;
}

////////////////////////////////////////////////////////////////////////////////
// Laplace3dBoundaryOperatorBase

template <typename Impl, typename BasisFunctionType, typename ResultType>
Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>::
Laplace3dBoundaryOperatorBase(
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label) :
    Base(domain, range, dualToRange, label), m_impl(new Impl),
    m_id(boost::make_shared<Laplace3dBoundaryOperatorId<BasisFunctionType> >(
             *this))
{
}

template <typename Impl, typename BasisFunctionType, typename ResultType>
Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>::
Laplace3dBoundaryOperatorBase(
        const Laplace3dBoundaryOperatorBase& other) :
    Base(other), m_impl(new Impl(*other.m_impl)), m_id(other.m_id)
{
}

template <typename Impl, typename BasisFunctionType, typename ResultType>
Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>::
~Laplace3dBoundaryOperatorBase()
{
}

template <typename Impl, typename BasisFunctionType, typename ResultType>
shared_ptr<const AbstractBoundaryOperatorId>
Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>::id() const
{
    return m_id;
}

template <typename Impl, typename BasisFunctionType, typename ResultType>
const typename Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>::
CollectionOfKernels&
Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>::
kernels() const
{
    return m_impl->kernels;
}

template <typename Impl, typename BasisFunctionType, typename ResultType>
const typename Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>::
CollectionOfBasisTransformations&
Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>::
testTransformations() const
{
    return m_impl->transformations;
}

template <typename Impl, typename BasisFunctionType, typename ResultType>
const typename Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>::
CollectionOfBasisTransformations&
Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>::
trialTransformations() const
{
    return m_impl->transformations;
}

template <typename Impl, typename BasisFunctionType, typename ResultType>
const typename Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>::
TestKernelTrialIntegral&
Laplace3dBoundaryOperatorBase<Impl, BasisFunctionType, ResultType>::
integral() const
{
    return m_impl->integral;
}

} // namespace Bempp

#endif
