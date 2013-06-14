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

#ifndef bempp_maxwell_3d_identity_operator_hpp
#define bempp_maxwell_3d_identity_operator_hpp

#include "../common/common.hpp"

#include "elementary_local_operator.hpp"

namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForOperators;
/** \endcond */

} // namespace Fiber

namespace Bempp
{

/** \cond FORWARD_DECL */
template <typename BasisFunctionType, typename ResultType> class BoundaryOperator;
/** \endcond */

/** \ingroup maxwell_3d
    \ingroup local_operators
    \brief "Identity operator" for Maxwell equations in 3D.

    This class represents an operator \f$I_{\mathrm{M}}\f$ whose weak form (under
    the standard inner product) is

    \f[
    \langle \boldsymbol u, I_{\mathrm{M}} \boldsymbol v \rangle =
    \int_\Gamma \boldsymbol u^*(\boldsymbol x) \cdot
    [\boldsymbol v(\boldsymbol x) \times \boldsymbol n(\boldsymbol x)] \,
    \mathrm{d}\Gamma(\boldsymbol x),
    \f]

    where \f$\boldsymbol n(\boldsymbol x)\f$ is the outward unit vector normal
    to the surface \f$\Gamma\f$ at the point \f$\boldsymbol x\f$. Thus, the
    weak form of this operator is identical to the weak form of the standard
    identity operator under the antisymmetric pseudo-inner product \f$\langle
    \cdot, \cdot\rangle_{\boldsymbol\tau, \Gamma}\f$ defined in \ref
    maxwell_3d.

    See AbstractBoundaryOperator for the documentation of the template
    parameters.

    Use the maxwell3dIdentityOperator() function to create a BoundaryOperator
    object wrapping a Maxwell3dIdentityOperator object.
 */
template <typename BasisFunctionType_, typename ResultType_>
class Maxwell3dIdentityOperator :
        public ElementaryLocalOperator<BasisFunctionType_, ResultType_>
{
    typedef ElementaryLocalOperator<BasisFunctionType_, ResultType_> Base;
public:
    /** \copydoc ElementaryLocalOperator::BasisFunctionType */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \copydoc ElementaryLocalOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc ElementaryLocalOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc ElementaryLocalOperator::QuadratureStrategy */
    typedef typename Base::QuadratureStrategy QuadratureStrategy;
    /** \copydoc ElementaryLocalOperator::LocalAssembler */
    typedef typename Base::LocalAssembler LocalAssembler;
    /** \copydoc ElementaryLocalOperator::CollectionOfBasisTransformations */
    typedef typename Base::CollectionOfBasisTransformations
    CollectionOfBasisTransformations;
    /** \copydoc ElementaryLocalOperator::TestTrialIntegral */
    typedef typename Base::TestTrialIntegral TestTrialIntegral;

    /** \brief Constructor.
     *
     *  \param[in] domain
     *    Function space being the domain of the operator.
     *  \param[in] range
     *    Function space being the range of the operator.
     *  \param[in] dualToRange
     *    Function space dual to the the range of the operator.
     *  \param[in] label
     *    Textual label of the operator. If empty, a unique label is generated
     *    automatically.
     *  \param[in] symmetry
     *    Symmetry of the weak form of the operator. Can be any combination of
     *    the flags defined in the enumeration type Symmetry.
     *
     *  All the three spaces must be defined on the same grid. */
    Maxwell3dIdentityOperator(const shared_ptr<const Space<BasisFunctionType> >& domain,
                     const shared_ptr<const Space<BasisFunctionType> >& range,
                     const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
                     const std::string& label = "",
                     int symmetry = NO_SYMMETRY);
    virtual ~Maxwell3dIdentityOperator();

private:
    virtual const CollectionOfBasisTransformations&
    testTransformations() const;

    virtual const CollectionOfBasisTransformations&
    trialTransformations() const;

    virtual const TestTrialIntegral& integral() const;

private:
    /** \cond PRIVATE */
    shared_ptr<TestTrialIntegral> m_integral;
    /** \endcond */
};

/** \relates Maxwell3dIdentityOperator
 *  \brief Construct a BoundaryOperator object wrapping a Maxwell3dIdentityOperator.
 *
 *  This convenience function constructs an abstract Maxwell identity operator
 *  and wraps it in a BoundaryOperator object.
 *
 *  \param[in] context
 *    A Context object that will be used to build the weak form of the
 *    identity operator when necessary.
 *  \param[in] domain
 *    Function space being the domain of the identity operator.
 *  \param[in] range
 *    Function space being the range of the identity operator.
 *  \param[in] dualToRange
 *    Function space dual to the the range of the identity operator.
 *  \param[in] label
 *    Textual label of the identity operator (optional, used for debugging).
 *  \param[in] symmetry
 *    Symmetry of the weak form of the operator. Can be any combination of
 *    the flags defined in the enumeration type Symmetry.
 *
 *  All the three spaces must be defined on the same grid. */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
maxwell3dIdentityOperator(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY);

} // namespace Bempp

#endif
