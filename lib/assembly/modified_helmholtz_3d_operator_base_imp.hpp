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

#ifndef bempp_modified_helmholtz_3d_operator_base_imp_hpp
#define bempp_modified_helmholtz_3d_operator_base_imp_hpp

#include "modified_helmholtz_3d_operator_base.hpp"

namespace Bempp
{

namespace
{

template <typename Impl>
inline typename Impl::KernelType waveNumberImpl(const Impl& impl)
{
    return impl.kernels.functor().waveNumber();
}

template <typename Impl>
inline void setWaveNumberImpl(Impl& impl, typename Impl::KernelType waveNumber)
{
    impl.kernels.functor().setWaveNumber(waveNumber);
}

} // namespace

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
ModifiedHelmholtz3dOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
ModifiedHelmholtz3dOperatorBase(
        const Space<BasisFunctionType>& domain,
        const Space<BasisFunctionType>& range,
        const Space<BasisFunctionType>& dualToRange,
        KernelType waveNumber,
        const std::string& label) :
    Base(domain, range, dualToRange, label), m_impl(new Impl(waveNumber))
{
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
ModifiedHelmholtz3dOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
ModifiedHelmholtz3dOperatorBase(
        const ModifiedHelmholtz3dOperatorBase& other) :
    Base(other), m_impl(new Impl(*other.m_impl))
{
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
ModifiedHelmholtz3dOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
~ModifiedHelmholtz3dOperatorBase()
{
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
typename ModifiedHelmholtz3dOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::KernelType
ModifiedHelmholtz3dOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
waveNumber() const
{
    return waveNumberImpl(*m_impl);
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
void
ModifiedHelmholtz3dOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
setWaveNumber(KernelType waveNumber)
{
    setWaveNumberImpl(*m_impl, waveNumber);
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
const typename ModifiedHelmholtz3dOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
CollectionOfKernels&
ModifiedHelmholtz3dOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
kernels() const
{
    return m_impl->kernels;
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
const typename ModifiedHelmholtz3dOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
CollectionOfBasisTransformations&
ModifiedHelmholtz3dOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
testTransformations() const
{
    return m_impl->transformations;
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
const typename ModifiedHelmholtz3dOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
CollectionOfBasisTransformations&
ModifiedHelmholtz3dOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
trialTransformations() const
{
    return m_impl->transformations;
}

template <typename Impl, typename BasisFunctionType,
          typename KernelType, typename ResultType>
const typename ModifiedHelmholtz3dOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
TestKernelTrialIntegral&
ModifiedHelmholtz3dOperatorBase<Impl, BasisFunctionType, KernelType, ResultType>::
integral() const
{
    return m_impl->integral;
}

} // namespace Bempp

#endif
