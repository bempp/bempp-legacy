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

#ifndef bempp_laplace_beltrami_3d_operator_hpp
#define bempp_laplace_beltrami_3d_operator_hpp

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

/** \relates LaplaceBeltrami3dOperator
 *  \brief Construct a BoundaryOperator object representing a Laplace-Beltrami
 *  operator in 3D.
 *
 *  The weak form of the Laplace-Beltrami operator \f$L\f$ in 3D is
 *
 *  \f[
 *  \langle \boldsymbol u, L \boldsymbol v \rangle =
 *  \int_\Gamma \boldsymbol \nabla_\Gamma u^*(\boldsymbol x) \cdot
 *  \boldsymbol \nabla_\Gamma v(\boldsymbol x)  \,
 *  \mathrm{d}\Gamma(\boldsymbol x),
 *  \f]
 *
 *  where \f$\boldsymbol \nabla_\Gamma\f$ is the surface gradient operator.
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
 *  All the three spaces must be defined on the same grid. See
 *  AbstractBoundaryOperator for the documentation of the template parameters.
 */
template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
laplaceBeltrami3dOperator(
        const shared_ptr<const Context<BasisFunctionType, ResultType> >& context,
        const shared_ptr<const Space<BasisFunctionType> >& domain,
        const shared_ptr<const Space<BasisFunctionType> >& range,
        const shared_ptr<const Space<BasisFunctionType> >& dualToRange,
        const std::string& label = "",
        int symmetry = NO_SYMMETRY);

} // namespace Bempp

#endif
