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

#ifndef bempp_modified_helmholtz_3d_hypersingular_boundaryoperator_hpp
#define bempp_modified_helmholtz_3d_hypersingular_boundaryoperator_hpp

#include "boundary_operator.hpp"
#include "helmholtz_3d_operators_common.hpp"
#include "../common/types.hpp"
#include "symmetry.hpp"

namespace Bempp {

/** \ingroup modified_helmholtz_3d
 *  \brief Construct a BoundaryOperator object representing the hypersingular
 *  operator associated with the modified Helmholtz equation in 3D.
 *
 *  \param[in] context
 *    A Context object that will be used to build the weak form of the
 *    boundary operator when necessary.
 *  \param[in] domain
 *    Function space being the domain of the boundary operator.
 *  \param[in] range
 *    Function space being the range of the boundary operator.
 *  \param[in] dualToRange
 *    Function space dual to the the range of the boundary operator.
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
 *  thrown.
 *
 *  If local-mode ACA assembly is requested (see AcaOptions::mode), after
 *  discretization, the weak form of this operator is stored as the product
 *
 *  \f[
 *     A = \sum_{i=1}^3 P_i A_{\textrm{d}} Q_i +
 *         k^2 \sum_{i=1}^3 R_i A_{\textrm{d}} S_i,
 *  \f]
 *
 *  where:
 *
 *  - \f$A_{\textrm{d}}\f$ is the weak form of the single-layer modified
 *    Helmholtz boundary operator discretised with test and trial functions
 *    being the restrictions of the basis functions of \p domain and \p range to
 *    individual elements;
 *
 *  - \f$P_i\f$ is the sparse matrix whose transpose represents the expansion
 *    of the \f$i\f$th component of the surface curl of the basis functions of
 *    \p dualToRange in the single-element test functions mentioned above;
 *
 *  - \f$Q_i\f$ is the sparse matrix representing the expansion of the \f$i\f$th
 *    component of the surface curl of the basis functions of \p domain in the
 *    single-element trial functions;
 *
 *  - \f$k\f$ is the wave number
 *
 *  - \f$R_i\f$ is the sparse matrix whose transpose represents the expansion
 *    of the basis functions of \p dualToRange, multiplied by the \f$i\f$th
 *    component of the local vector normal to the surface on which functions
 *    from \p dualToRange live, in the single-element test functions;
 *
 *  - \f$S_i\f$ is the sparse matrix representing the expansion of basis
 *    functions of \p domain, multiplied by the \f$i\f$th component of the local
 *    vector normal to the surface on which functions from \p domain live, in
 *    the single-element trial functions.
 *
 *  \tparam BasisFunctionType
 *    Type of the values of the basis functions into
 *    which functions acted upon by the operator are expanded.
 *  \tparam KernelType
 *    Type of the values of the kernel functions occurring
 *    in the integrand of the operator.
 *  \tparam ResultType
 *    Type used to represent elements of the weak form form of the operator.
 *
 *  The latter three template parameters can take the following values: \c
 *  float, \c double, <tt>std::complex<float></tt> and
 *  <tt>std::complex<double></tt>. All types must have the same precision: for
 *  instance, mixing \c float with <tt>std::complex<double></tt> is not
 *  allowed. The parameter \p ResultType is by default set to "larger" of \p
 *  BasisFunctionType and \p KernelType, e.g. for \p BasisFunctionType = \c
 *  double and \p KernelType = <tt>std::complex<double></tt> it is set to
 *  <tt>std::complex<double></tt>. You should override that only if you set
 *  both \p BasisFunctionType and \p KernelType to a real type, but you want
 *  the entries of the operator's weak form to be stored as complex numbers.
 *
 *  Note that setting \p KernelType to a real type implies that the wave number
 *  must also be chosen purely real.
 */
template <typename BasisFunctionType, typename KernelType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
modifiedHelmholtz3dHypersingularBoundaryOperator(
    const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    KernelType waveNumber, const std::string &label = "",
    int symmetry = NO_SYMMETRY, bool useInterpolation = false,
    int interpPtsPerWavelength = DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY,
    const BoundaryOperator<BasisFunctionType, ResultType> &externalSlp =
        BoundaryOperator<BasisFunctionType, ResultType>());


template <typename BasisFunctionType, typename KernelType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
modifiedHelmholtz3dHypersingularBoundaryOperator(
    const ParameterList& parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    KernelType waveNumber, const std::string &label = "",
    int symmetry = NO_SYMMETRY, bool useInterpolation = false,
    int interpPtsPerWavelength = DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY,
    const BoundaryOperator<BasisFunctionType, ResultType> &externalSlp =
        BoundaryOperator<BasisFunctionType, ResultType>());

} // namespace Bempp

#endif
