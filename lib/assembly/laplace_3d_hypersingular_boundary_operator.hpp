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

#ifndef bempp_laplace_3d_hypersingular_boundary_operator_hpp
#define bempp_laplace_3d_hypersingular_boundary_operator_hpp

#include "boundary_operator.hpp"
#include "symmetry.hpp"

namespace Bempp
{

/** \ingroup laplace_3d
 *  \brief Construct a BoundaryOperator object representing the hypersingular
 *  boundary operator associated with the Laplace equation in 3D.
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
 *  thrown. */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplace3dHypersingularBoundaryOperator(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY);

/** \ingroup laplace_3d
 *  \brief Construct a "synthetic" representation of the the hypersingular
 *  boundary operator associated with the Laplace equation in 3D.
 *
 *  This function creates a hypersingular Laplace boundary operator \f$\mathcal
 *  A\f$ whose weak form is stored as
 *
 *  \[
 *     A = \sum_{\alpha=1}^3 P_\alpha A_{\textrm{d}} Q_\alpha,
 *  \]
 *
 *  where \f$A_{\textrm{d}}\f$ is the weak form of the single-layer Laplace
 *  boundary operator discretised with the <em>discontinuous</em> test and
 *  trial function spaces passed in the parameters \p internalTestSpace and \p
 *  internalTrialSpace, \f$Q_alpha\f$ is the sparse matrix representing the
 *  \f$\alpha\f$th component of the surface curl of the basis functions of \p
 *  domain in the basis of \p internalTrialSpace, and \f$P\f$ is the sparse
 *  matrix whose transpose represents the \f$\alpha\f$th component of the
 *  surface curl of the basis functions of \p dualTorange in the basis of \p
 *  internalTestSpace.
 *
 *  The discrete weak form of this operator is usually assembled faster than
 *  that of the "standard" operator created by
 *  laplace3dHypersingularBoundaryOperator() (ACA is faster for operators
 *  discretised with discontinuous spaces), but uses up more memory. In
 *  addition, multiplication of the discrete weak form of this operator by a
 *  vector is approximately three times slower than for the discrete weak form
 *  of the "standard" operator.
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
 *  an exception is thrown. */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplace3dSyntheticHypersingularBoundaryOperator(
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
