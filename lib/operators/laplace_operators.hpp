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

#ifndef bempp_laplace_operators_hpp
#define bempp_laplace_operators_hpp

#include "../assembly/symmetry.hpp"
#include "../common/shared_ptr.hpp"
#include "../common/types.hpp"

namespace Bempp {

template <typename BasisFunctionType, typename KernelType, typename ResultType>
class ElementaryIntegralOperator;
template <typename BasisFunctionType, typename ResultType>
class PotentialOperator;
template <typename BasisFunctionType> class Space;

template <typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<
    const ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>>
laplaceSingleLayerBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label = "", int symmetry = NO_SYMMETRY);

template <typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<
    const ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>>
laplaceDoubleLayerBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label = "", int symmetry = NO_SYMMETRY);

template <typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<
    const ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>>
laplaceAdjointDoubleLayerBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label = "", int symmetry = NO_SYMMETRY);

template <typename BasisFunctionType, typename KernelType, typename ResultType>
shared_ptr<
    const ElementaryIntegralOperator<BasisFunctionType, KernelType, ResultType>>
laplaceHypersingularBoundaryOperator(
    const ParameterList &parameterList,
    const shared_ptr<const Space<BasisFunctionType>> &domain,
    const shared_ptr<const Space<BasisFunctionType>> &range,
    const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
    const std::string &label = "", int symmetry = NO_SYMMETRY);

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>>
laplaceSingleLayerPotentialOperator(
    const shared_ptr<const Space<BasisFunctionType>> &space,
    const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>
        &evaluationPoints,
    const ParameterList &parameterList);

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>>
laplaceDoubleLayerPotentialOperator(
    const shared_ptr<const Space<BasisFunctionType>> &space,
    const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>
        &evaluationPoints,
    const ParameterList &parameterList);

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>>
laplaceSingleLayerGradientPotentialOperator(
    const shared_ptr<const Space<BasisFunctionType>> &space,
    const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>
        &evaluationPoints,
    const ParameterList &parameterList);

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>>
laplaceDoubleLayerGradientPotentialOperator(
    const shared_ptr<const Space<BasisFunctionType>> &space,
    const Matrix<typename Fiber::ScalarTraits<BasisFunctionType>::RealType>
        &evaluationPoints,
    const ParameterList &parameterList);
}

#include "laplace_operators_imp.hpp"

#endif
