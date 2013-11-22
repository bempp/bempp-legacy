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

#include "fmm_black_box_double_layer.hpp"
#include "fmm_black_box.hpp"

#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

template <typename KernelType, typename ValueType>
void FmmBlackBoxDoubleLayer<KernelType, ValueType>::evaluateTrial(
            const arma::Col<CoordinateType>& point,
            const arma::Col<CoordinateType>& normal,
            const arma::Col<CoordinateType>& multipole, // [-1,1]
            const arma::Col<CoordinateType>& nodeCentre,
            const arma::Col<CoordinateType>& nodeSize,
            arma::Col<ValueType>& result) const
{
    this->evaluateAtGaussPointDiffS(point, normal, multipole, 
        nodeCentre, nodeSize, result);
}

template <typename KernelType, typename ValueType>
void FmmBlackBoxDoubleLayer<KernelType, ValueType>::evaluateTest(
            const arma::Col<CoordinateType>& point,
            const arma::Col<CoordinateType>& normal,
            const arma::Col<CoordinateType>& multipole,
            const arma::Col<CoordinateType>& nodeCentre,
            const arma::Col<CoordinateType>& nodeSize,
            arma::Col<ValueType>& result) const
{
    this->evaluateAtGaussPointS(point, normal, multipole, 
        nodeCentre, nodeSize, result);
}

// should be templated on KernelType and ResultType, but not added to explicit 
// instantiation yet. The following is equivalent.
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(FmmBlackBoxDoubleLayer);

} // namespace Bempp
