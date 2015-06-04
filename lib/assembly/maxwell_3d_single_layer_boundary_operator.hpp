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

#ifndef bempp_maxwell_3d_single_layer_boundary_operator_hpp
#define bempp_maxwell_3d_single_layer_boundary_operator_hpp

#include "boundary_operator.hpp"
#include "helmholtz_3d_operators_common.hpp"
#include "symmetry.hpp"
#include "../common/types.hpp"

#include "../common/scalar_traits.hpp"

namespace Bempp {

/** \ingroup maxwell_3d
 *  \brief Construct a single-layer boundary operator associated with the
 *  Maxwell equations in 3D.
 *
 *  This function constructs a BoundaryOperator object representing
 *  the single-layer boundary operator \f$\boldsymbol S\f$ for the
 *  Maxwell equations in 3D, as defined in \ref maxwell_3d.
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
 *    Wave number. See \ref maxwell_3d for its definition.
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
 *  This operator supports local-mode ACA assembly. It does not support
 *  hybrid-mode assembly, though.
 *
 *  If local-mode ACA assembly is requested (see AcaOptions::mode), after
 *  discretization, the weak form of this operator is stored as
 *
 *  \f[
 *     A = -ik \sum_{i=1}^3 P_i A_{\textrm{d}} Q_i -
 *         \frac{1}{ik} R_i A_{\textrm{d}} S_i,
 *  \f]
 *
 *  where:
 *
 *  - \f$k\f$ is the wave number
 *
 *  - \f$A_{\textrm{d}}\f$ is the weak form of the single-layer modified
 *    Helmholtz boundary operator discretised with test and trial functions
 *    being the restrictions of the basis functions of \p domain and \p range to
 *    individual elements;
 *
 *  - \f$P_i\f$ is the sparse matrix whose transpose represents the expansion of
 *    the \f$i\f$th component of the basis functions of \p dualToRange in the
 *    single-element test functions mentioned above;
 *
 *  - \f$Q_i\f$ is the sparse matrix representing the expansion of the \f$i\f$th
 *    component of the basis functions of \p domain in the single-element trial
 *    functions;
 *
 *  - \f$R_i\f$ is the sparse matrix whose transpose represents the expansion of
 *    the surface divergence of the basis functions of \p dualToRange in the
 *    single-element test functions;
 *
 *  - \f$S_i\f$ is the sparse matrix representing the expansion of the surface
 *    divergence of the basis functions of \p domain in the single-element trial
 *    functions.
 *
 *  \tparam BasisFunctionType
 *    Type of the values of the basis functions into which functions acted upon
 *    by the operator are expanded. It can take the following values: \c float,
 *    \c double, <tt>std::complex<float></tt> and
 *    <tt>std::complex<double></tt>.
 *
 *  \see maxwell3dSingleLayerPotentialOperator()
 */
template <typename BasisFunctionType>
BoundaryOperator<BasisFunctionType,
                 typename ScalarTraits<BasisFunctionType>::ComplexType>
maxwell3dSingleLayerBoundaryOperator(
    const shared_ptr<const Context<
        BasisFunctionType,
        typename ScalarTraits<BasisFunctionType>::ComplexType>> &context,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    typename ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
    const std::string &label = "", int symmetry = NO_SYMMETRY,
    bool useInterpolation = false,
    int interpPtsPerWavelength = DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY);

template <typename BasisFunctionType>
BoundaryOperator<BasisFunctionType,
                 typename ScalarTraits<BasisFunctionType>::ComplexType>
maxwell3dSingleLayerBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    typename ScalarTraits<BasisFunctionType>::ComplexType waveNumber,
    const std::string &label = "", int symmetry = NO_SYMMETRY,
    bool useInterpolation = false,
    int interpPtsPerWavelength = DEFAULT_HELMHOLTZ_INTERPOLATION_DENSITY);
} // namespace Bempp

#endif
