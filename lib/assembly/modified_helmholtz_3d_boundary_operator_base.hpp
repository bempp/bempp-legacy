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

#ifndef bempp_modified_helmholtz_3d_boundary_operator_base_hpp
#define bempp_modified_helmholtz_3d_boundary_operator_base_hpp

#include "elementary_singular_integral_operator.hpp"

#include <boost/scoped_ptr.hpp>

#include "helmholtz_3d_operators_common.hpp"

namespace Bempp
{

/** \cond FORWARD_DECL */
template <typename BasisFunctionType> class ModifiedHelmholtz3dBoundaryOperatorId;
/** \endcond */

/** \ingroup helmholtz_3d
 *  \brief Base class for boundary operators for the modified Helmholtz
 *  equation in 3D.
 *
 *  \tparam Impl
 *    Type of the internal implementation object.
 *  \tparam BasisFunctionType_
 *    Type of the values of the basis functions into
 *    which functions acted upon by the operator are expanded.
 *  \tparam KernelType_
 *    Type of the values of the kernel functions occurring
 *    in the integrand of the operator.
 *  \tparam ResultType_
 *    Type used to represent elements of the weak form form of the operator.
 *
 *  The latter three template parameters can take the following values: \c
 *  float, \c double, <tt>std::complex<float></tt> and
 *  <tt>std::complex<double></tt>. All types must have the same precision: for
 *  instance, mixing \c float with <tt>std::complex<double></tt> is not
 *  allowed. The parameter \p ResultType_ is by default set to "larger" of \p
 *  BasisFunctionType_ and \p KernelType_, e.g. for \p BasisFunctionType_ = \c
 *  double and \p KernelType_ = <tt>std::complex<double></tt> it is set to
 *  <tt>std::complex<double></tt>. You should override that only if you set
 *  both \p BasisFunctionType_ and \p KernelType_ to a real type, but you want
 *  the entries of the operator's weak form to be stored as complex numbers.
 *
 *  Note that setting \p KernelType_ to a real type implies that the wave number
 *  must also be chosen purely real.
 *
 *  \see modified_helmholtz_3d
 */
template <typename Impl, typename BasisFunctionType_,
          typename KernelType_, typename ResultType_>
class ModifiedHelmholtz3dBoundaryOperatorBase :
        public ElementarySingularIntegralOperator<
        BasisFunctionType_, KernelType_, ResultType_>
{
    typedef ElementarySingularIntegralOperator<
    BasisFunctionType_, KernelType_, ResultType_>
    Base;

    friend class ModifiedHelmholtz3dBoundaryOperatorId<BasisFunctionType_>;

public:
    /** \brief Type of the values of the basis functions into which functions
     *  acted upon by the operator are expanded. */
    typedef typename Base::BasisFunctionType BasisFunctionType;
    /** \brief Type of the values of kernel functions. */
    typedef typename Base::KernelType KernelType;
    /** \copydoc ElementaryIntegralOperator::ResultType */
    typedef typename Base::ResultType ResultType;
    /** \copydoc ElementaryIntegralOperator::CoordinateType */
    typedef typename Base::CoordinateType CoordinateType;
    /** \copydoc ElementaryIntegralOperator::CollectionOfBasisTransformations */
    typedef typename Base::CollectionOfBasisTransformations
    CollectionOfBasisTransformations;
    /** \copydoc ElementaryIntegralOperator::CollectionOfKernels */
    typedef typename Base::CollectionOfKernels CollectionOfKernels;
    /** \copydoc ElementaryIntegralOperator::TestKernelTrialIntegral */
    typedef typename Base::TestKernelTrialIntegral TestKernelTrialIntegral;

    /** \brief Constructor.
     *
     *  \param[in] domain
     *    Function space being the domain of the operator.
     *  \param[in] range
     *    Function space being the range of the operator.
     *  \param[in] dualToRange
     *    Function space dual to the the range of the operator.
     *  \param[in] waveNumber
     *    Wave number. See \ref modified_helmholtz_3d for its definition.
     *  \param[in] label
     *    Textual label of the operator. If empty, a unique label is generated
     *    automatically.
     *  \param[in] symmetry
     *    Symmetry of the weak form of the operator. Can be any combination of the
     *    flags defined in the enumeration type Symmetry.
     *  \param[in] useInterpolation
     *    If set to \p false (default), the standard exp() function will be used to
     *    evaluate the exponential factor occurring in the kernel. If set to \p
     *    true, the exponential factor will be evaluated by piecewise-cubic
     *    interpolation of values calculated in advance on a regular grid. This
     *    normally speeds up calculations, but might result in a loss of accuracy.
     *    This is an experimental feature: use it at your own risk.
     *  \param[in] interpPtsPerWavelength
     *    If \p useInterpolation is set to \p true, this parameter determines the
     *    number of points per "effective wavelength" (defined as \f$2\pi/|k|\f$,
     *    where \f$k\f$ = \p waveNumber) used to construct the interpolation grid.
     *    The default value (5000) is normally enough to reduce the relative or
     *    absolute error, *whichever is smaller*, below 100 * machine precision. If
     *    \p useInterpolation is set to \p false, this parameter is ignored.
     *
     *  None of the shared pointers may be null and the spaces \p range and \p
     *  dualToRange must be defined on the same grid, otherwise an exception is
     *  thrown. */
    ModifiedHelmholtz3dBoundaryOperatorBase(
            const shared_ptr<const Space<BasisFunctionType> >& domain,
            const shared_ptr<const Space<BasisFunctionType> >& range,
            const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
            KernelType waveNumber,
            const std::string& label = "",
            int symmetry = NO_SYMMETRY,
            bool useInterpolation = false,
            int interpPtsPerWavelength = DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY);
    /** \brief Copy constructor. */
    ModifiedHelmholtz3dBoundaryOperatorBase(
            const ModifiedHelmholtz3dBoundaryOperatorBase& other);

    /** \brief Destructor. */
    virtual ~ModifiedHelmholtz3dBoundaryOperatorBase();

    /** \brief Return the wave number set previously in the constructor. */
    KernelType waveNumber() const;

    /** \brief Return the identifier of this operator.
     *
     *  Two boundary operators related to the modified Helmholtz equation are
     *  treated as identical, and hence having the same weak form, if they have
     *  the same C++ type (e.g.
     *  <tt>Helmholtz3dDoubleLayerBoundaryOperator<double, double, double></tt>),
     *  domain space, range space, space dual to range and wave number. */
    virtual shared_ptr<const AbstractBoundaryOperatorId> id() const;

private:
    virtual const CollectionOfKernels& kernels() const;
    virtual const CollectionOfBasisTransformations&
    testTransformations() const;
    virtual const CollectionOfBasisTransformations&
    trialTransformations() const;
    virtual const TestKernelTrialIntegral& integral() const;

private:
    /** \cond PRIVATE */
    boost::scoped_ptr<Impl> m_impl;
    shared_ptr<AbstractBoundaryOperatorId> m_id;
    /** \endcond */
};

} // namespace Bempp

#endif
