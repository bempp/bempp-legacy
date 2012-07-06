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
Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
Helmholtz3dBoundaryOperatorBase(       
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        KernelType waveNumber,
        const std::string& label) :
    Base(domain, range, dualToRange, label), m_impl(new Impl(waveNumber))
{
}

template <typename Impl, typename BasisFunctionType>
Helmholtz3dBoundaryOperatorBase<Impl, BasisFunctionType>::
Helmholtz3dBoundaryOperatorBase(
        const Helmholtz3dBoundaryOperatorBase& other) :
    Base(other), m_impl(new Impl(*other.m_impl))
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
