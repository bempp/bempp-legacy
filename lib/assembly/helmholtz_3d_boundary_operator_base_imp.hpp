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

#ifndef bempp_helmholtz_3d_boundary_operator_base_imp_hpp
#define bempp_helmholtz_3d_boundary_operator_base_imp_hpp

#include "helmholtz_3d_boundary_operator_base.hpp"
#include "abstract_boundary_operator_id.hpp"
#include "../common/boost_make_shared_fwd.hpp"

namespace Bempp
{

namespace
{

template <typename Impl>
inline typename Impl::KernelType waveNumberImpl(const Impl& impl)
{
    return impl.kernels.functor().waveNumber() *
            typename Impl::KernelType(0., 1.);
}

template <typename Impl>
inline void setWaveNumberImpl(Impl& impl, typename Impl::KernelType waveNumber)
{
    impl.kernels.functor().setWaveNumber(
                waveNumber / typename Impl::KernelType(0., 1.));
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
// Helmholtz3dBoundaryOperatorId

template <typename BasisFunctionType>
class Helmholtz3dBoundaryOperatorId : public AbstractBoundaryOperatorId
{
public:
    template <typename Impl>
    explicit Helmholtz3dBoundaryOperatorId(
            const Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>& op);
    virtual size_t hash() const;
    virtual void dump() const;
    virtual bool isEqual(const AbstractBoundaryOperatorId &other) const;

private:
    const std::type_info& m_typeInfo;
    const Space<BasisFunctionType>* m_domain;
    const Space<BasisFunctionType>* m_range;
    const Space<BasisFunctionType>* m_dualToRange;
    typename ScalarTraits<BasisFunctionType>::ComplexType m_waveNumber;
};

template <typename BasisFunctionType>
template <typename Impl>
Helmholtz3dBoundaryOperatorId<BasisFunctionType>::Helmholtz3dBoundaryOperatorId(
        const Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>& op) :
    m_typeInfo(typeid(op)),
    m_domain(op.domain().get()), m_range(op.range().get()),
    m_dualToRange(op.dualToRange().get()),
    m_waveNumber(op.waveNumber())
{
}

template <typename BasisFunctionType>
size_t Helmholtz3dBoundaryOperatorId<BasisFunctionType>::hash() const
{
    size_t result = tbb::tbb_hasher(m_typeInfo.name());
    tbb_hash_combine(result, m_domain);
    tbb_hash_combine(result, m_range);
    tbb_hash_combine(result, m_dualToRange);
    tbb_hash_combine(result, std::abs(m_waveNumber));
    return result;
}

template <typename BasisFunctionType>
void Helmholtz3dBoundaryOperatorId<BasisFunctionType>::dump() const
{
    std::cout << m_typeInfo.name() << ", " << m_domain << ", "
              << m_range << ", " << m_dualToRange 
              << ", " << m_waveNumber << std::endl;
}

template <typename BasisFunctionType>
bool Helmholtz3dBoundaryOperatorId<BasisFunctionType>::isEqual(
        const AbstractBoundaryOperatorId &other) const
{
    // dynamic_cast won't suffice since we want to make sure both objects
    // are of exactly the same type (dynamic_cast would succeed for a subclass)
    if (typeid(other) == typeid(*this)) {
        const Helmholtz3dBoundaryOperatorId& otherCompatible =
            static_cast<const Helmholtz3dBoundaryOperatorId&>(other);
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
// Helmholtz3dBoundaryOperatorBase

template <typename Impl, typename BasisFunctionType>
Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
Helmholtz3dBoundaryOperatorBase(       
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        KernelType waveNumber,
        const std::string& label) :
    Base(domain, range, dualToRange, label), m_impl(new Impl(waveNumber)),
    m_id(boost::make_shared<Helmholtz3dBoundaryOperatorId<BasisFunctionType> >(
             *this))
{
}

template <typename Impl, typename BasisFunctionType>
Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
Helmholtz3dBoundaryOperatorBase(
        const Helmholtz3dBoundaryOperatorBase& other) :
    Base(other), m_impl(new Impl(*other.m_impl)), m_id(other.m_id)
{
}

template <typename Impl, typename BasisFunctionType>
Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
~Helmholtz3dBoundaryOperatorBase()
{
}

template <typename Impl, typename BasisFunctionType>
typename Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::KernelType
Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
waveNumber() const
{
    return waveNumberImpl(*m_impl);
}

template <typename Impl, typename BasisFunctionType>
void
Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
setWaveNumber(KernelType waveNumber)
{
    setWaveNumberImpl(*m_impl, waveNumber);
}

template <typename Impl, typename BasisFunctionType>
shared_ptr<const AbstractBoundaryOperatorId>
Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::id() const
{
    return m_id;
}

template <typename Impl, typename BasisFunctionType>
const typename Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
CollectionOfKernels&
Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
kernels() const
{
    return m_impl->kernels;
}

template <typename Impl, typename BasisFunctionType>
const typename Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
CollectionOfBasisTransformations&
Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
testTransformations() const
{
    return m_impl->transformations;
}

template <typename Impl, typename BasisFunctionType>
const typename Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
CollectionOfBasisTransformations&
Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
trialTransformations() const
{
    return m_impl->transformations;
}

template <typename Impl, typename BasisFunctionType>
const typename Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
TestKernelTrialIntegral&
Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
integral() const
{
    return m_impl->integral;
}

} // namespace Bempp

#endif
