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

#ifndef bempp_helmholtz_3d_potential_operator_base_hpp
#define bempp_helmholtz_3d_potential_operator_base_hpp

#include "elementary_potential_operator.hpp"

#include <boost/scoped_ptr.hpp>

namespace Bempp
{

/** \ingroup helmholtz_3d
 *  \brief Base class for boundary operators for the Helmholtz equation in 3D.
 *
 *  \tparam Impl
 *    Type of the internal implementation pointer.
 *  \tparam BasisFunctionType
 *    Type used to represent the values of basis functions. It can take the
 *    following values: \c float, \c double, <tt>std::complex<float></tt> and
 *    <tt>std::complex<double></tt>.
 *
 *  \see helmholtz_3d */
template <typename Impl, typename BasisFunctionType_>
class Helmholtz3dPotentialOperatorBase :
        public ElementaryPotentialOperator<
        BasisFunctionType_,
        typename ScalarTraits<BasisFunctionType_>::ComplexType,
        typename ScalarTraits<BasisFunctionType_>::ComplexType>
{
    typedef ElementaryPotentialOperator<
    BasisFunctionType_,
    typename ScalarTraits<BasisFunctionType_>::ComplexType,
    typename ScalarTraits<BasisFunctionType_>::ComplexType>
    Base;
public:
    typedef typename Base::BasisFunctionType BasisFunctionType;
    typedef typename Base::KernelType KernelType;
    typedef typename Base::ResultType ResultType;
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::CollectionOfBasisTransformations
    CollectionOfBasisTransformations;
    typedef typename Base::CollectionOfKernels CollectionOfKernels;
    typedef typename Base::KernelTrialIntegral KernelTrialIntegral;

    Helmholtz3dPotentialOperatorBase(KernelType waveNumber);
    Helmholtz3dPotentialOperatorBase(const Helmholtz3dPotentialOperatorBase& other);
    virtual ~Helmholtz3dPotentialOperatorBase();

    void setWaveNumber(KernelType waveNumber);
    KernelType waveNumber() const;

private:
    virtual const CollectionOfKernels& kernels() const;
    virtual const CollectionOfBasisTransformations&
    trialTransformations() const;
    virtual const KernelTrialIntegral& integral() const;

private:
    boost::scoped_ptr<Impl> m_impl;
};

} // namespace Bempp

#endif
