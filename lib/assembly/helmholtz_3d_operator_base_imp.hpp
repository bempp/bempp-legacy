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

#ifndef bempp_helmholtz_3d_operator_base_imp_hpp
#define bempp_helmholtz_3d_operator_base_imp_hpp

#include "helmholtz_3d_operator_base.hpp"

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

template <typename Impl, typename BasisFunctionType>
Helmholtz3dOperatorBase<Impl, BasisFunctionType>::
Helmholtz3dOperatorBase(       
        const Space<BasisFunctionType>& domain,
        const Space<BasisFunctionType>& range,
        const Space<BasisFunctionType>& dualToRange,
        KernelType waveNumber,
        const std::string& label) :
    Base(domain, range, dualToRange, label), m_impl(new Impl(waveNumber))
{
}

template <typename Impl, typename BasisFunctionType>
Helmholtz3dOperatorBase<Impl, BasisFunctionType>::
Helmholtz3dOperatorBase(
        const Helmholtz3dOperatorBase& other) :
    Base(other), m_impl(new Impl(*other.m_impl))
{
}

template <typename Impl, typename BasisFunctionType>
Helmholtz3dOperatorBase<Impl, BasisFunctionType>::
~Helmholtz3dOperatorBase()
{
}

template <typename Impl, typename BasisFunctionType>
typename Helmholtz3dOperatorBase<Impl, BasisFunctionType>::KernelType
Helmholtz3dOperatorBase<Impl, BasisFunctionType>::
waveNumber() const
{
    return waveNumberImpl(*m_impl);
}

template <typename Impl, typename BasisFunctionType>
void
Helmholtz3dOperatorBase<Impl, BasisFunctionType>::
setWaveNumber(KernelType waveNumber)
{
    setWaveNumberImpl(*m_impl, waveNumber);
}

template <typename Impl, typename BasisFunctionType>
const typename Helmholtz3dOperatorBase<Impl, BasisFunctionType>::
CollectionOfKernels&
Helmholtz3dOperatorBase<Impl, BasisFunctionType>::
kernels() const
{
    return m_impl->kernels;
}

template <typename Impl, typename BasisFunctionType>
const typename Helmholtz3dOperatorBase<Impl, BasisFunctionType>::
CollectionOfBasisTransformations&
Helmholtz3dOperatorBase<Impl, BasisFunctionType>::
testTransformations() const
{
    return m_impl->transformations;
}

template <typename Impl, typename BasisFunctionType>
const typename Helmholtz3dOperatorBase<Impl, BasisFunctionType>::
CollectionOfBasisTransformations&
Helmholtz3dOperatorBase<Impl, BasisFunctionType>::
trialTransformations() const
{
    return m_impl->transformations;
}

template <typename Impl, typename BasisFunctionType>
const typename Helmholtz3dOperatorBase<Impl, BasisFunctionType>::
TestKernelTrialIntegral&
Helmholtz3dOperatorBase<Impl, BasisFunctionType>::
integral() const
{
    return m_impl->integral;
}

} // namespace Bempp

#endif
