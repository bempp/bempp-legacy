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

#ifndef bempp_maxwell_3d_single_layer_potential_operator_hpp
#define bempp_maxwell_3d_single_layer_potential_operator_hpp

#include "helmholtz_3d_potential_operator_base.hpp"

namespace Bempp
{

/** \cond PRIVATE */
template <typename BasisFunctionType>
struct Maxwell3dSingleLayerPotentialOperatorImpl;
/** \endcond */

/** \ingroup maxwell_3d
 *  \brief Single-layer potential operator for the Maxwell's equations in 3D.
 *
    This class represents the single-layer potential operator occurring in the
    Stratton-Chu theorem (Buffa and Hiptmair, Theorem 6) on the representation
    of solutions of the time-harmonic Maxwell equations reduced to the form

    \f[
    \boldsymbol\nabla \times \boldsymbol\nabla \times \boldsymbol u -
    k^2 \boldsymbol u = 0,
    \f],

    where \f$u\f$ stands either for the electric or the magnetic field, in a
    homogeneous region. The theorem states that any solution \f$u\f$ of the
    above equation in a bounded domain \f$\Omega\f$ with boundary \f$\Gamma\f$
    has the representation

    \f[ \boldsymbol u =
    \boldsymbol \Psi_{\mathrm{SL}}(\gamma_{\mathrm{N,int}}\boldsymbol u) +
    \boldsymbol \Psi_{\mathrm{DL}}(\gamma_{\mathrm{D,int}}\boldsymbol u).
    \f]

    (all the symbols will be defined shortly). Furthermore, any solution \f$u\f$ satisfying Maxwell's equations in a domain \f$\Omega^{\mathrm c}\f$ extending to infinity as well as the Silver-Mueller radiation conditions at infinity has the representation

    \f[ \boldsymbol u =
    -\boldsymbol \Psi_{\mathrm{SL}}(\gamma_{\mathrm{N,ext}}\boldsymbol u) -
    \boldsymbol \Psi_{\mathrm{DL}}(\gamma_{\mathrm{D,ext}}\boldsymbol u).
    \f]

    The action of the <em>single-layer potential operator</em> is defined by

    \f[
    (\boldsymbol \Psi_{\mathrm{SL}} \boldsymbol v)(\boldsymbol x) \equiv
    \mathrm i k \int_{\Gamma} G(\boldsymbol x, \boldsymbol y)
    \boldsymbol v(\boldsymbol y) \mathrm\Gamma(\boldsymbol y) -
    \frac{1}{\mathrm i k}
    \boldsymbol\nabla_{\boldsymbol x}
    \int_{\Gamma} G(\boldsymbol x, \boldsymbol y)
    (\boldsymbol \nabla_\Gamma \cdot \boldsymbol v)(\boldsymbol y)
    \mathrm\Gamma(\boldsymbol y),
    \f]

    where

    \f[
    G(\boldsymbol x, \boldsymbol y) \equiv
    \frac{\exp(\mathrm i k \lvert \boldsymbol x - \boldsymbol y \rvert)}
    {4\pi \lvert \boldsymbol x - \boldsymbol y \rvert}
    \f]

    is the Green's function of the Helmholtz equation. The action of the
    <em>double-layer potential operator</em>, in turn, is defined by
    \f[
    (\Psi_{\mathrm{SL}} \boldsymbol v)(\boldsymbol x) \equiv
    \boldsymbol\nabla_{\boldsymbol x} \times
    \int_{\Gamma} G(\boldsymbol x, \boldsymbol y)
    \boldsymbol v(\boldsymbol y) \mathrm\Gamma(\boldsymbol y).
    \f]

    The <em>interior Dirichlet trace</em> \f$\gamma_{\mathrm{D,int}}\boldsymbol
    u\f$ at a point \f$\boldsymbol x \in \Gamma\f$ is defined as

    \f[(\gamma_{\mathrm{D,int}}\boldsymbol u)(\boldsymbol x) \equiv
    \boldsymbol u\rvert_{\Gamma,\mathrm{int}}(\boldsymbol x) \times
    \boldsymbol n(\boldsymbol x),\f]

    where \f$\boldsymbol n\f$
    is the outward unit vector normal to \f$\Gamma\f$ at \f$\boldsymbol x\f$ and
    \f$\boldsymbol u\rvert_{\Gamma,\mathrm{int}}(\boldsymbol x)\f$
    is the limit of \f$\boldsymbol u(\boldsymbol y)\f$ as
    \f$\boldsymbol y\f$ approaches \f$\boldsymbol x\f$
    from within \f$\Omega\f$.
    The <em>interior Neumann trace</em> \f$\gamma_{\mathrm{N,int}}\boldsymbol
    u\f$ at \f$\boldsymbol x \in \Gamma\f$ is defined as

    \f[(\gamma_{\mathrm{N,int}}\boldsymbol u)(\boldsymbol x) \equiv
    \frac{1}{\mathrm i k}
    (\boldsymbol \nabla \times \boldsymbol u)\rvert_{\Gamma,\mathrm{int}}
    (\boldsymbol x) \times
    \boldsymbol n(\boldsymbol x).\f]

    Note that compared to Buffa and Hiptmair, we include an additional
    \f$\mathrm i\f$ factor in the denominator of the Neumann trace, as this
    makes both the Neumann trace and the potential operators real-valued in the
    special case of \f$k\f$ purely imaginary (and \f$\boldsymbol u\f$ purely real).

    The <em>exterior</em> Dirichlet and Neumann traces are defined analogously,
    with the relevant quantities assumed to approach the point \f$\boldsymbol x\f$
    from within the complement of \f$\Omega\f$.


    \tparam BasisFunctionType_
      Type of the values of the basis functions into which functions acted upon
      by the operator are expanded. It can take the following values: \c float,
      \c double, <tt>std::complex<float></tt> and
      <tt>std::complex<double></tt>.

    \see maxwell_3d */
template <typename BasisFunctionType_>
class Maxwell3dSingleLayerPotentialOperator :
        public Helmholtz3dPotentialOperatorBase<
        Maxwell3dSingleLayerPotentialOperatorImpl<BasisFunctionType_>,
        BasisFunctionType_>
{
    typedef Helmholtz3dPotentialOperatorBase<
    Maxwell3dSingleLayerPotentialOperatorImpl<BasisFunctionType_>,
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

    /** \brief Constructor.
     *
     *  \param[in] waveNumber
     *    Wave number. See \ref maxwell_3d for its definition. */
    Maxwell3dSingleLayerPotentialOperator(KernelType waveNumber);
    /** \copydoc Helmholtz3dPotentialOperatorBase::~Helmholtz3dPotentialOperatorBase */
    virtual ~Maxwell3dSingleLayerPotentialOperator();
};

} // namespace Bempp

#endif
