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

#include "fmm_modified_helmholtz_3d_double_layer_high_frequency.hpp"
#include "fmm_modified_helmholtz_3d_high_frequency.hpp"

#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

template <typename ValueType>
void FmmModifiedHelmholtz3dDoubleLayerHighFrequency<ValueType>::evaluateTrial(
            const arma::Col<CoordinateType>& point,
            const arma::Col<CoordinateType>& normal,
            const arma::Col<CoordinateType>& khat,
            const arma::Col<CoordinateType>& nodeCentre,
            const arma::Col<CoordinateType>& nodeSize,
            arma::Col<ValueType>& result) const
{
    ValueType kappa = FmmModifiedHelmholtz3dHighFrequency<ValueType>::kappa();
    arma::Col<CoordinateType> r = nodeCentre - point;
    result(0) =  kappa*exp( -kappa*dot(khat, r) )*dot(khat, normal);
}

template <typename ValueType>
void FmmModifiedHelmholtz3dDoubleLayerHighFrequency<ValueType>::evaluateTest(
            const arma::Col<CoordinateType>& point,
            const arma::Col<CoordinateType>& normal,
            const arma::Col<CoordinateType>& khat,
            const arma::Col<CoordinateType>& nodeCentre,
            const arma::Col<CoordinateType>& nodeSize,
            arma::Col<ValueType>& result) const
{
    ValueType kappa = FmmModifiedHelmholtz3dHighFrequency<ValueType>::kappa();
    arma::Col<CoordinateType> r = point - nodeCentre;
    result(0) =  exp( -kappa*dot(khat, r) );
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(FmmModifiedHelmholtz3dDoubleLayerHighFrequency);

} // namespace Bempp
