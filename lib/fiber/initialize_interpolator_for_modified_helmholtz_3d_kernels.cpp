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

#include "initialize_interpolator_for_modified_helmholtz_3d_kernels.hpp"
#include "explicit_instantiation.hpp"

namespace Fiber
{

template <typename ValueType>
void initializeInterpolatorForModifiedHelmholtz3dKernels(
        ValueType waveNumber,
        typename ScalarTraits<ValueType>::RealType maxDist,
        int interpPtsPerWavelength,
        HermiteInterpolator<ValueType>& interpolator)
{
    typedef typename ScalarTraits<ValueType>::RealType CoordinateType;
    const CoordinateType minDist = 0.;
    const CoordinateType wavelength = 2. * M_PI / std::abs(waveNumber);
    const int pointCount =
            (maxDist - minDist) / wavelength * interpPtsPerWavelength + 1;
    std::vector<ValueType> values(pointCount), derivatives(pointCount);
    for (int i = 0; i < pointCount; ++i) {
        CoordinateType dist =
            minDist + (maxDist - minDist) * i / CoordinateType(pointCount - 1);
        ValueType exponential = exp(-waveNumber * dist);
        values[i] = exponential;
        derivatives[i] = -waveNumber * exponential;
    }
    interpolator.initialize(minDist, maxDist, values, derivatives);
}

#define INSTANTIATE_FUNCTION(KERNEL) \
   template void initializeInterpolatorForModifiedHelmholtz3dKernels( \
       KERNEL, ScalarTraits<KERNEL>::RealType, int, \
       HermiteInterpolator<KERNEL>&);

FIBER_ITERATE_OVER_KERNEL_TYPES(INSTANTIATE_FUNCTION);

} // namespace Fiber
