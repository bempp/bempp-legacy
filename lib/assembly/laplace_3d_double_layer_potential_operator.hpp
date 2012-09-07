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

#ifndef bempp_laplace_3d_double_layer_potential_operator_hpp
#define bempp_laplace_3d_double_layer_potential_operator_hpp

#include "laplace_3d_potential_operator_base.hpp"

namespace Bempp
{

/** \cond PRIVATE */
template <typename BasisFunctionType, typename ResultType>
struct Laplace3dDoubleLayerPotentialOperatorImpl;
/** \endcond */

/** \ingroup laplace_3d
 *  \brief Double-layer potential operator for the Laplace equation in 3D.
 *
 *  \tparam BasisFunctionType_
 *    Type of the values of the basis functions into
 *    which functions acted upon by the operator are expanded.
 *  \tparam ResultType_
 *    Type of the values of the potential.
 *
 *  Both template parameters can take the following values: \c float, \c
 *  double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>.
 *  Both types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. The parameter \p ResultType_
 *  is by default set to \p BasisFunctionType_. You should override that only if
 *  you set \p BasisFunctionType_ to a real type, but you want the values of the
 *  potential to be stored as complex numbers.
 *
 *  \see laplace_3d */
template <typename BasisFunctionType_, typename ResultType_ = BasisFunctionType_>
class Laplace3dDoubleLayerPotentialOperator:
        public Laplace3dPotentialOperatorBase<
        Laplace3dDoubleLayerPotentialOperatorImpl<BasisFunctionType_, ResultType_>,
        BasisFunctionType_,
        ResultType_>
{
    typedef Laplace3dPotentialOperatorBase<
    Laplace3dDoubleLayerPotentialOperatorImpl<BasisFunctionType_, ResultType_>,
    BasisFunctionType_,
    ResultType_>
    Base;
public:
    /** \copydoc Laplace3dPotentialOperatorBase::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc Laplace3dPotentialOperatorBase::KernelType */
    typedef typename Base::KernelType KernelType;
    /** \copydoc Laplace3dPotentialOperatorBase::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc Laplace3dPotentialOperatorBase::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc Laplace3dPotentialOperatorBase::CollectionOfBasisTransformations */
    typedef typename Base::CollectionOfBasisTransformations
    CollectionOfBasisTransformations;
    /** \copydoc Laplace3dPotentialOperatorBase::CollectionOfKernels */
    typedef typename Base::CollectionOfKernels CollectionOfKernels;
    /** \copydoc Laplace3dPotentialOperatorBase::KernelTrialIntegral */
    typedef typename Base::KernelTrialIntegral KernelTrialIntegral;

    /** \copydoc Laplace3dPotentialOperatorBase::Laplace3dPotentialOperatorBase */
    Laplace3dDoubleLayerPotentialOperator();
    /** \copydoc Laplace3dPotentialOperatorBase::~Laplace3dPotentialOperatorBase */
    virtual ~Laplace3dDoubleLayerPotentialOperator();
};

} // namespace Bempp

#endif
