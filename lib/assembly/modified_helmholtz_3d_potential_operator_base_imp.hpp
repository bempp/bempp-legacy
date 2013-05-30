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

#ifndef bempp_modified_helmholtz_3d_potential_operator_base_imp_hpp
#define bempp_modified_helmholtz_3d_potential_operator_base_imp_hpp

#include "modified_helmholtz_3d_potential_operator_base.hpp"

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

template <typename Impl, typename BasisFunctionType>
ModifiedHelmholtz3dPotentialOperatorBase<Impl, BasisFunctionType>::
ModifiedHelmholtz3dPotentialOperatorBase(KernelType waveNumber) :
    m_impl(new Impl(waveNumber))
{
}

template <typename Impl, typename BasisFunctionType>
ModifiedHelmholtz3dPotentialOperatorBase<Impl, BasisFunctionType>::
ModifiedHelmholtz3dPotentialOperatorBase(const ModifiedHelmholtz3dPotentialOperatorBase& other) :
    m_impl(new Impl(*other.m_impl))
{
}

template <typename Impl, typename BasisFunctionType>
ModifiedHelmholtz3dPotentialOperatorBase<Impl, BasisFunctionType>::
~ModifiedHelmholtz3dPotentialOperatorBase()
{
}

template <typename Impl, typename BasisFunctionType>
ModifiedHelmholtz3dPotentialOperatorBase<Impl, BasisFunctionType>&
ModifiedHelmholtz3dPotentialOperatorBase<Impl, BasisFunctionType>::
operator=(const ModifiedHelmholtz3dPotentialOperatorBase& rhs)
{
    if (this != &rhs) {
        Base::operator=(rhs);
        m_impl.reset(new Impl(*rhs.m_impl));
    }
}

template <typename Impl, typename BasisFunctionType>
typename ModifiedHelmholtz3dPotentialOperatorBase<Impl, BasisFunctionType>::KernelType
ModifiedHelmholtz3dPotentialOperatorBase<Impl, BasisFunctionType>::
waveNumber() const
{
    return waveNumberImpl(*m_impl);
}

template <typename Impl, typename BasisFunctionType>
const typename ModifiedHelmholtz3dPotentialOperatorBase<Impl, BasisFunctionType>::
CollectionOfKernels&
ModifiedHelmholtz3dPotentialOperatorBase<Impl, BasisFunctionType>::
kernels() const
{
    return m_impl->kernels;
}

template <typename Impl, typename BasisFunctionType>
const typename ModifiedHelmholtz3dPotentialOperatorBase<Impl, BasisFunctionType>::
CollectionOfBasisTransformations&
ModifiedHelmholtz3dPotentialOperatorBase<Impl, BasisFunctionType>::
trialTransformations() const
{
    return m_impl->transformations;
}

template <typename Impl, typename BasisFunctionType>
const typename ModifiedHelmholtz3dPotentialOperatorBase<Impl, BasisFunctionType>::
KernelTrialIntegral&
ModifiedHelmholtz3dPotentialOperatorBase<Impl, BasisFunctionType>::
integral() const
{
    return m_impl->integral;
}

} // namespace Bempp

#endif
