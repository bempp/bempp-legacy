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

#ifndef bempp_laplace_3d_single_layer_boundary_operator_hpp
#define bempp_laplace_3d_single_layer_boundary_operator_hpp

#include "boundary_operator.hpp"
#include "symmetry.hpp"

namespace Bempp
{

/** \ingroup laplace_3d
 *  \brief Construct a BoundaryOperator object representing the
 *  single-layer boundary operator associated with the Laplace equation in 3D.
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
 *  \param[in] label
 *    Textual label of the operator. If empty, a unique label is generated
 *    automatically.
 *  \param[in] symmetry
 *    Symmetry of the weak form of the operator. Can be any combination of the
 *    flags defined in the enumeration type Symmetry.
 *
 *  None of the shared pointers may be null and the spaces \p range and \p
 *  dualToRange must be defined on the same grid, otherwise an exception is
 *  thrown.
 *
 *  \tparam BasisFunctionType
 *    Type of the values of the basis functions into
 *    which functions acted upon by the operator are expanded.
 *  \tparam ResultType
 *    Type used to represent elements of the weak form of the operator.
 *
 *  Both template parameters can take the following values: \c float, \c
 *  double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>.
 *  Both types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. The parameter \p ResultType
 *  is by default set to \p BasisFunctionType. You should override that only if
 *  you set \p BasisFunctionType to a real type, but you want the entries of
 *  the operator's weak form to be stored as complex numbers.
 */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplace3dSingleLayerBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY);

/** \ingroup laplace_3d
 *  \brief Construct a "synthetic" representation of the single-layer boundary
 *  operator associated with the Laplace equation in 3D.
 *
 *  This function creates a single-layer Laplace boundary operator \f$\mathcal
 *  A\f$ whose weak form is stored as the product
 *
 *  \f[
 *     A = P A_{\textrm{d}} Q,
 *  \f]
 *
 *  where \f$A_{\textrm{d}}\f$ is the weak form of \f$\mathcal A\f$ discretised
 *  with the <em>discontinuous</em> test and trial function spaces passed in the
 *  parameters \p internalTestSpace and \p internalTrialSpace, \f$Q\f$ is the
 *  sparse matrix mapping the degrees of freedom of \p domain to those of \p
 *  internalTrialSpace, and \f$P\f$ the matrix mapping the degrees of freedom of
 *  \p internalTestSpace to \p dualToRange.
 *
 *  If <tt>internalTrialSpace == domain</tt> or <tt>internalTestSpace ==
 *  dualToRange</tt>, \f$P\f$ or \f$Q\f$ becomes an identity matrix and is
 *  not physically created.
 *
 *  In ACA mode, synthetic operators are, as a rule, built more quickly than
 *  "standard" operators (ACA is faster for operators discretised with
 *  discontinuous spaces), but use up more memory and their matrix-vector
 *  product may be more time-consuming.
 *
 *  See the documentation of SyntheticIntegralOperator for more information
 *  about the concept of synthetic operators.
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
 *  \param[in] internalTrialSpace
 *    Trial function space used in the discretisation of \f$\mathcal A\f$ to
 *    \f$A_{\textrm{d}}\f$. It must be a discontinuous space, with basis
 *    functions extending over single elements only.
 *  \param[in] internalTestSpace
 *    Test function space used in the discretisation of \f$\mathcal A\f$ to
 *    \f$A_{\textrm{d}}\f$. It must be a discontinuous space, with basis
 *    functions extending over single elements only.
 *  \param[in] label
 *    Textual label of the operator. If empty, a unique label is generated
 *    automatically.
 *  \param[in] symmetry
 *    Symmetry of the weak form of the operator. Can be any combination of the
 *    flags defined in the enumeration type Symmetry.
 *
 *  None of the shared pointers may be null. \p internalTestSpace must be
 *  defined on the same grid as \p range and \p dualToRange, while \p
 *  internalTrialSpace must be defined on the same grid as \p domain; otherwise
 *  an exception is thrown.
 *
 *   *  \tparam BasisFunctionType
 *    Type of the values of the basis functions into
 *    which functions acted upon by the operator are expanded.
 *  \tparam ResultType
 *    Type used to represent elements of the weak form of the operator.
 *
 *  Both template parameters can take the following values: \c float, \c
 *  double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>.
 *  Both types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. The parameter \p ResultType
 *  is by default set to \p BasisFunctionType. You should override that only if
 *  you set \p BasisFunctionType to a real type, but you want the entries of
 *  the operator's weak form to be stored as complex numbers.
 */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplace3dSyntheticSingleLayerBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const shared_ptr<const Space<BasisFunctionType> >& internalTrialSpace,
        const shared_ptr<const Space<BasisFunctionType> >& internalTestSpace,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY);

} // namespace Bempp

#endif
