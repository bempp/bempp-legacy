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
#ifndef laplace_3d_calderon_projector_hpp
#define laplace_3d_calderon_projector_hpp

#include "blocked_boundary_operator.hpp"
#include "../fiber/scalar_traits.hpp"
#include "../common/types.hpp"
#include "helmholtz_3d_operators_common.hpp"
#include "symmetry.hpp"

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
BlockedBoundaryOperator<BasisFunctionType, ResultType>
laplace3dExteriorCalderonProjector(
    const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
    const shared_ptr<const Space<BasisFunctionType>> &hminusSpace,
    const shared_ptr<const Space<BasisFunctionType>> &hplusSpace,
    const std::string &label = "");

template <typename BasisFunctionType, typename ResultType>
BlockedBoundaryOperator<BasisFunctionType, ResultType>
laplace3dInteriorCalderonProjector(
    const shared_ptr<const Context<BasisFunctionType, ResultType>> &context,
    const shared_ptr<const Space<BasisFunctionType>> &hminusSpace,
    const shared_ptr<const Space<BasisFunctionType>> &hplusSpace,
    const std::string &label = "");

template <typename BasisFunctionType, typename ResultType>
BlockedBoundaryOperator<BasisFunctionType, ResultType>
laplace3dExteriorCalderonProjector(
    const ParameterList& parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &hminusSpace,
    const shared_ptr<const Space<BasisFunctionType>> &hplusSpace,
    const std::string &label = "");

template <typename BasisFunctionType, typename ResultType>
BlockedBoundaryOperator<BasisFunctionType, ResultType>
laplace3dInteriorCalderonProjector(
    const ParameterList& parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &hminusSpace,
    const shared_ptr<const Space<BasisFunctionType>> &hplusSpace,
    const std::string &label = "");


} // namespace Bempp

#endif // laplace_3d_calderon_projector_hpp
