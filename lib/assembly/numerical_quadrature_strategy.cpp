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

#include "numerical_quadrature_strategy.hpp"

#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/numerical_quadrature_strategy_imp.hpp"

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
NumericalQuadratureStrategy<BasisFunctionType,
                            ResultType>::NumericalQuadratureStrategy()
    : Base() {}

template <typename BasisFunctionType, typename ResultType>
NumericalQuadratureStrategy<BasisFunctionType, ResultType>::
    NumericalQuadratureStrategy(const AccuracyOptionsEx &accuracyOptions)
    : Base(accuracyOptions) {}

template <typename BasisFunctionType, typename ResultType>
NumericalQuadratureStrategy<BasisFunctionType, ResultType>::
    NumericalQuadratureStrategy(
        const AccuracyOptionsEx &accuracyOptions,
        const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>>
            &singleQuadratureRuleFamily,
        const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType>>
            &doubleQuadratureRuleFamily)
    : Base(accuracyOptions, singleQuadratureRuleFamily,
           doubleQuadratureRuleFamily) {}

template <typename BasisFunctionType, typename ResultType>
NumericalQuadratureStrategy<BasisFunctionType, ResultType>::
    NumericalQuadratureStrategy(
        const shared_ptr<const QuadratureDescriptorSelectorFactory<
            BasisFunctionType>> &quadratureDescriptorSelectorFactory,
        const shared_ptr<const SingleQuadratureRuleFamily<CoordinateType>>
            &singleQuadratureRuleFamily,
        const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType>>
            &doubleQuadratureRuleFamily)
    : Base(quadratureDescriptorSelectorFactory, singleQuadratureRuleFamily,
           doubleQuadratureRuleFamily) {}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(
    NumericalQuadratureStrategy);

} // namespace Bempp
