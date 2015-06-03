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
#include "../common/global_parameters.hpp"

namespace Bempp {

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
 *  If local-mode ACA assembly is requested (see AcaOptions::mode), after
 *  discretization, the weak form of this operator is stored as the product
 *
 *  \f[
 *     P A_{\textrm{d}} Q,
 *  \f]
 *
 *  where \f$A_{\textrm{d}}\f$ is the weak form of this operator discretized
 *  with test and trial functions being the restrictions of the basis functions
 *  of \p domain and \p range to individual elements; \f$Q\f$ is the sparse
 *  matrix representing the expansion of the basis functions of \p domain in the
 *  just mentioned single-element trial functions; and \f$P\f$ is the sparse
 *  matrix whose transpose represents the expansion of the basis functions of \p
 *  dualToRange in the single-element test functions.
 *
 *  If each basis function of \p domain extends over a single element only
 *  (i.e. <tt>domain->isDiscontinuous()</tt> method returns \c true), the matrix
 *  \f$Q\f$ is omitted from the discretized weak form. The same applies to \p
 *  dualToDomain and the matrix \f$P\f$.
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
    const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label = "", int symmetry = NO_SYMMETRY);

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplace3dSingleLayerBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label = "", int symmetry = NO_SYMMETRY);
} // namespace Bempp

#endif
