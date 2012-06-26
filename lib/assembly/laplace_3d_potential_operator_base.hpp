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

#ifndef bempp_laplace_3d_potential_operator_base_hpp
#define bempp_laplace_3d_potential_operator_base_hpp

#include "elementary_potential_operator.hpp"

#include <boost/scoped_ptr.hpp>

namespace Bempp
{

/** \ingroup laplace_3d
 *  \brief Base class for potentials related to the Laplace equation in 3D.
 *
 *  \tparam Impl
 *    Type of the internal implementation pointer.
 *  \tparam BasisFunctionType
 *    Type used to represent the values of basis functions.
 *  \tparam ResultType
 *    Type used to represent entries in the discrete form of the operator.
 *
 *  The latter two template parameters can take the following values: \c float, \c
 *  double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>.
 *  Both types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. The parameter \p ResultType
 *  is by default set to \p BasisFunctionType. You should override that only if
 *  you set \p BasisFunctionType to a real type, but you want the potential
 *  values to be stored as complex numbers.
 *
 *  \see laplace_3d */
template <typename Impl,
          typename BasisFunctionType_, typename ResultType_ = BasisFunctionType_>
class Laplace3dPotentialOperatorBase :
        public ElementaryPotentialOperator<
        BasisFunctionType_,
        typename ScalarTraits<ResultType_>::RealType,
        ResultType_>
{
    typedef ElementaryPotentialOperator<
    BasisFunctionType_,
    typename ScalarTraits<ResultType_>::RealType,
    ResultType_>
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

    Laplace3dPotentialOperatorBase();
    Laplace3dPotentialOperatorBase(const Laplace3dPotentialOperatorBase& other);
    virtual ~Laplace3dPotentialOperatorBase();

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
