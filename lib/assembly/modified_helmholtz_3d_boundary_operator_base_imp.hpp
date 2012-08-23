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

#ifndef bempp_modified_helmholtz_3d_boundary_operator_base_imp_hpp
#define bempp_modified_helmholtz_3d_boundary_operator_base_imp_hpp

#include "modified_helmholtz_3d_boundary_operator_base.hpp"
#include "abstract_boundary_operator_id.hpp"
#include "../common/boost_make_shared_fwd.hpp"

namespace Bempp
{

namespace
{

template <typename Impl>
inline typename Impl::KernelType waveNumberImpl(const Impl& impl)
{
    return impl.kernels.functor().waveNumber();
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
// ModifiedHelmholtz3dBoundaryOperatorId

template <typename BasisFunctionType>
class ModifiedHelmholtz3dBoundaryOperatorId : public AbstractBoundaryOperatorId
{
public:
    template <typename Impl, typename KernelType, typename ResultType>
    explicit ModifiedHelmholtz3dBoundaryOperatorId(
            const ModifiedHelmholtz3dBoundaryOperatorBase<
            Impl, BasisFunctionType, KernelType, ResultType>& op);
    virtual size_t hash() const;
    virtual void dump() const;
    virtual bool isEqual(const AbstractBoundaryOperatorId &other) const;

private:
    const std::type_info& m_typeInfo;
    const Space<BasisFunctionType>* m_domain;
    const Space<BasisFunctionType>* m_range;
    const Space<BasisFunctionType>* m_dualToRange;
    // it is simplest to just always store the wave number as a complex number
    typename ScalarTraits<BasisFunctionType>::ComplexType m_waveNumber;
};

template <typename BasisFunctionType>
template <typename Impl, typename KernelType, typename ResultType>
ModifiedHelmholtz3dBoundaryOperatorId<BasisFunctionType>::
ModifiedHelmholtz3dBoundaryOperatorId(
        const ModifiedHelmholtz3dBoundaryOperatorBase<
        Impl, BasisFunctionType, KernelType, ResultType>& op) :
    m_typeInfo(typeid(op)),
    m_domain(op.domain().get()), m_range(op.range().get()),
    m_dualToRange(op.dualToRange().get()),
    m_waveNumber(op.waveNumber())
{
}

template <typename BasisFunctionType>
size_t ModifiedHelmholtz3dBoundaryOperatorId<BasisFunctionType>::hash() const
{
    size_t result = tbb::tbb_hasher(m_typeInfo.name());
    tbb_hash_combine(result, m_domain);
    tbb_hash_combine(result, m_range);
    tbb_hash_combine(result, m_dualToRange);
    tbb_hash_combine(result, std::abs(m_waveNumber));
    return result;
}

template <typename BasisFunctionType>
void ModifiedHelmholtz3dBoundaryOperatorId<BasisFunctionType>::dump() const
{
    std::cout << m_typeInfo.name() << ", " << m_domain << ", "
              << m_range << ", " << m_dualToRange 
              << ", " << m_waveNumber << std::endl;
}

template <typename BasisFunctionType>
bool ModifiedHelmholtz3dBoundaryOperatorId<BasisFunctionType>::isEqual(
        const AbstractBoundaryOperatorId &other) const
{
    // dynamic_cast won't suffice since we want to make sure both objects
    // are of exactly the same type (dynamic_cast would succeed for a subclass)
    if (typeid(other) == typeid(*this)) {
        const ModifiedHelmholtz3dBoundaryOperatorId& otherCompatible =
            static_cast<const ModifiedHelmholtz3dBoundaryOperatorId&>(other);
        return (m_typeInfo == otherCompatible.m_typeInfo &&
                m_domain == otherCompatible.m_domain &&
                m_range == otherCompatible.m_range &&
                m_dualToRange == otherCompatible.m_dualToRange &&
                m_waveNumber == otherCompatible.m_waveNumber);
    }
    else
        return false;
}

////////////////////////////////////////////////////////////////////////////////
// ModifiedHelmholtz3dBoundaryOperatorBase

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
ModifiedHelmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
ModifiedHelmholtz3dBoundaryOperatorBase(
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        KernelType waveNumber,
        const std::string& label) :
    Base(domain, range, dualToRange, label), m_impl(new Impl(waveNumber)),
    m_id(boost::make_shared<ModifiedHelmholtz3dBoundaryOperatorId<BasisFunctionType> >(
             *this))
{
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
ModifiedHelmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
ModifiedHelmholtz3dBoundaryOperatorBase(
        const ModifiedHelmholtz3dBoundaryOperatorBase& other) :
    Base(other), m_impl(new Impl(*other.m_impl)), m_id(other.m_id)
{
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
ModifiedHelmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
~ModifiedHelmholtz3dBoundaryOperatorBase()
{
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
typename ModifiedHelmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::KernelType
ModifiedHelmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
waveNumber() const
{
    return waveNumberImpl(*m_impl);
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
shared_ptr<const AbstractBoundaryOperatorId>
ModifiedHelmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
id() const
{
    return m_id;
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
const typename ModifiedHelmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
CollectionOfKernels&
ModifiedHelmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
kernels() const
{
    return m_impl->kernels;
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
const typename ModifiedHelmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
CollectionOfBasisTransformations&
ModifiedHelmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
testTransformations() const
{
    return m_impl->transformations;
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
const typename ModifiedHelmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
CollectionOfBasisTransformations&
ModifiedHelmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
trialTransformations() const
{
    return m_impl->transformations;
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
const typename ModifiedHelmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
TestKernelTrialIntegral&
ModifiedHelmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
integral() const
{
    return m_impl->integral;
}

} // namespace Bempp

#endif
