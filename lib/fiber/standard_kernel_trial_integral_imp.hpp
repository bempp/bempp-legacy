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

#ifndef fiber_standard_standard_kernel_trial_integral_imp_hpp
#define fiber_standard_standard_kernel_trial_integral_imp_hpp

#include "standard_kernel_trial_integral.hpp"

#include "geometrical_data.hpp"

#include <algorithm>

namespace Fiber
{

/*
template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class IntegrandFunctor
{
public:
    typedef BasisFunctionType_ BasisFunctionType;
    typedef KernelType_ KernelType;
    typedef ResultType_ ResultType;
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    void addGeometricalDependencies(int& trialGeomDeps) const;

    ResultType evaluate(
            const ConstGeometricalDataSlice<CoordinateType>& trialGeomData,
            const CollectionOf2dSlicesOfConst4dArrays<KernelType>& kernels,
            const CollectionOf1dSlicesOfConst2dArrays<BasisFunctionType>&
            weightedTrialTransformations) const;
};
*/

template <typename IntegrandFunctor>
void StandardKernelTrialIntegral<IntegrandFunctor>::
addGeometricalDependencies(size_t& trialGeomDeps) const
{
    m_functor.addGeometricalDependencies(trialGeomDeps);
}

template <typename IntegrandFunctor>
int StandardKernelTrialIntegral<IntegrandFunctor>::resultDimension() const
{
    return m_functor.resultDimension();
}

template <typename IntegrandFunctor>
void StandardKernelTrialIntegral<IntegrandFunctor>::evaluate(
        const GeometricalData<CoordinateType>& trialGeomData,
        const CollectionOf4dArrays<KernelType>& kernels,
        const CollectionOf2dArrays<ResultType>& weightedTrialTransformations,
        _2dArray<ResultType>& result) const
{
    const size_t evalPointCount = kernels[0].extent(2);
    const size_t quadPointCount = kernels[0].extent(3);
    const int resultComponentCount = m_functor.resultDimension();

    for (size_t i = 0; i < kernels.size(); ++i) {
        assert(kernels[i].extent(2) == evalPointCount);
        assert(kernels[i].extent(3) == quadPointCount);
    }
    for (size_t i = 0; i < weightedTrialTransformations.size(); ++i) {
        assert((int)weightedTrialTransformations[i].extent(0) == resultComponentCount);
        assert(weightedTrialTransformations[i].extent(1) == quadPointCount);
    }

    result.set_size(resultComponentCount, evalPointCount);
    std::fill(result.begin(), result.end(), 0.);

    for (size_t evalPoint = 0; evalPoint < evalPointCount; ++evalPoint)
        for (int dim = 0; dim < resultComponentCount; ++dim)
            for (size_t quadPoint = 0; quadPoint < quadPointCount; ++quadPoint)
                result(dim, evalPoint) += m_functor.evaluate(
                            trialGeomData.const_slice(quadPoint),
                            kernels.const_slice(evalPoint, quadPoint),
                            weightedTrialTransformations.const_slice(quadPoint));
}

} // namespace Fiber

#endif
