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

#ifndef bempp_helmholtz_3d_far_field_single_layer_potential_operator_hpp
#define bempp_helmholtz_3d_far_field_single_layer_potential_operator_hpp

#include "helmholtz_3d_potential_operator_base.hpp"

namespace Bempp
{

/** \cond PRIVATE */
template <typename BasisFunctionType>
struct Helmholtz3dFarFieldSingleLayerPotentialOperatorImpl;
/** \endcond */

/** \ingroup helmholtz_3d
 *  \brief Part of the far-field pattern of a radiating solution of the
 *  Helmholtz equation in 3D.
 *
 *  The far-field pattern of a radiating solution to the Helmholtz equation in
 *  3D is given by (Colton and Kress, p. 21)
 *  \f[
 *      u_\infty(\hat x) =
 *          -u_{\infty,\mathrm{SL}}(\hat x) + u_{\infty,\mathrm{DL}}(\hat x),
 *  \f]
 *  where \f$\hat x\f$ is a vector lying on the unit sphere,
 *  \f[
 *      u_{\infty,\mathrm{SL}} = \frac{1}{4\pi} \int_S
 *          e^{-i k \hat x \cdot y}
 *          \frac{\partial u}{\partial n(y)}\,
 *          dS(y)
 *  \f]
 *  and
 *  \f[
 *      u_{\infty,\mathrm{DL}} = \frac{1}{4\pi} \int_S
 *          \frac{\partial e^{-i k \hat x \cdot y}}{\partial n_(y)}
 *          u(y)\,
 *          dS(y),
 *  \f]
 *  with \f$n(y)\f$ denoting the unit vector normal to \e S at \e y.
 *
 *  This operator, applied to a grid function \f$v(y) = \frac{\partial
 *  u}{\partial n_y}(y) \f$ defined on a surface \e S, and evaluated at a point
 *  \f$\hat x\f$ lying on the unit sphere, yields
 *  \f$u_{\infty,\mathrm{SL}}(\hat x)\f$ as defined above.
 *
 *  \tparam BasisFunctionType_
 *    Type of the values of the basis functions into which functions acted upon
 *    by the operator are expanded. It can take the following values: \c float,
 *    \c double, <tt>std::complex<float></tt> and
 *    <tt>std::complex<double></tt>.
 *
 *  \see helmholtz_3d */
template <typename BasisFunctionType_>
class Helmholtz3dFarFieldSingleLayerPotentialOperator :
        public Helmholtz3dPotentialOperatorBase<
        Helmholtz3dFarFieldSingleLayerPotentialOperatorImpl<BasisFunctionType_>,
        BasisFunctionType_>
{
    typedef Helmholtz3dPotentialOperatorBase<
    Helmholtz3dFarFieldSingleLayerPotentialOperatorImpl<BasisFunctionType_>,
    BasisFunctionType_>
    Base;
public:
    /** \copydoc Helmholtz3dPotentialOperatorBase::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc Helmholtz3dPotentialOperatorBase::KernelType */
    typedef typename Base::KernelType KernelType;
    /** \copydoc Helmholtz3dPotentialOperatorBase::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc Helmholtz3dPotentialOperatorBase::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc Helmholtz3dPotentialOperatorBase::CollectionOfBasisTransformations */
    typedef typename Base::CollectionOfBasisTransformations
    CollectionOfBasisTransformations;
    /** \copydoc Helmholtz3dPotentialOperatorBase::CollectionOfKernels */
    typedef typename Base::CollectionOfKernels CollectionOfKernels;
    /** \copydoc Helmholtz3dPotentialOperatorBase::KernelTrialIntegral */
    typedef typename Base::KernelTrialIntegral KernelTrialIntegral;

    /** \copydoc Helmholtz3dPotentialOperatorBase::Helmholtz3dPotentialOperatorBase */
    Helmholtz3dFarFieldSingleLayerPotentialOperator(KernelType waveNumber);
    /** \copydoc Helmholtz3dPotentialOperatorBase::~Helmholtz3dPotentialOperatorBase */
    virtual ~Helmholtz3dFarFieldSingleLayerPotentialOperator();
};

} // namespace Bempp

#endif
