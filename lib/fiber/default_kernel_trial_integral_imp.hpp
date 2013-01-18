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

#ifndef fiber_default_kernel_trial_integral_imp_hpp
#define fiber_default_kernel_trial_integral_imp_hpp

#include "default_kernel_trial_integral.hpp"

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

    void evaluate(
            const ConstGeometricalDataSlice<CoordinateType>& trialGeomData,
            const CollectionOf2dSlicesOfConst4dArrays<KernelType>& kernels,
            const CollectionOf1dSlicesOfConst2dArrays<BasisFunctionType>&
            weightedTrialTransformations,
            std::vector<ResultType>& value) const;
};
*/

template <typename IntegrandFunctor>
void DefaultKernelTrialIntegral<IntegrandFunctor>::
addGeometricalDependencies(size_t& trialGeomDeps) const
{
    trialGeomDeps |= INTEGRATION_ELEMENTS;
    m_functor.addGeometricalDependencies(trialGeomDeps);
}

template <typename IntegrandFunctor>
int DefaultKernelTrialIntegral<IntegrandFunctor>::resultDimension() const
{
    return m_functor.resultDimension();
}

template <typename IntegrandFunctor>
void DefaultKernelTrialIntegral<IntegrandFunctor>::evaluate(
        const GeometricalData<CoordinateType>& trialGeomData,
        const CollectionOf4dArrays<KernelType>& kernels,
        const CollectionOf2dArrays<ResultType>& trialTransformations,
        const std::vector<CoordinateType>& weights,
        _2dArray<ResultType>& result) const
{
    const size_t evalPointCount = kernels[0].extent(2);
    const size_t quadPointCount = kernels[0].extent(3);
    const int resultComponentCount = m_functor.resultDimension();

    for (size_t i = 0; i < kernels.size(); ++i) {
        assert(kernels[i].extent(2) == evalPointCount);
        assert(kernels[i].extent(3) == quadPointCount);
    }
    for (size_t i = 0; i < trialTransformations.size(); ++i) {
        assert(trialTransformations[i].extent(1) == quadPointCount);
    }
    assert(weights.size() == quadPointCount);

    result.set_size(resultComponentCount, evalPointCount);
    std::fill(result.begin(), result.end(), 0.);

    std::vector<ResultType> value(resultComponentCount);

    for (size_t evalPoint = 0; evalPoint < evalPointCount; ++evalPoint)
        for (size_t quadPoint = 0; quadPoint < quadPointCount; ++quadPoint) {
            m_functor.evaluate(
                        trialGeomData.const_slice(quadPoint),
                        kernels.const_slice(evalPoint, quadPoint),
                        trialTransformations.const_slice(quadPoint),
                        value);
            for (int dim = 0; dim < resultComponentCount; ++dim)
                result(dim, evalPoint) += value[dim] * weights[quadPoint];
        }
}

template <typename IntegrandFunctor>
void DefaultKernelTrialIntegral<IntegrandFunctor>::evaluateWithPureWeights(
        const GeometricalData<CoordinateType>& trialGeomData,
        const CollectionOf4dArrays<KernelType>& kernels,
        const CollectionOf3dArrays<BasisFunctionType>& trialTransformations,
        const std::vector<CoordinateType>& weights,
        _3dArray<ResultType>& result) const
{
    const size_t trialDofCount = trialTransformations[0].extent(1);
    const size_t evalPointCount = kernels[0].extent(2);
    const size_t quadPointCount = kernels[0].extent(3);
    const int resultComponentCount = m_functor.resultDimension();

    for (size_t i = 0; i < kernels.size(); ++i) {
        assert(kernels[i].extent(2) == evalPointCount);
        assert(kernels[i].extent(3) == quadPointCount);
    }
    for (size_t i = 0; i < trialTransformations.size(); ++i) {
        assert(trialTransformations[i].extent(1) == trialDofCount);
        assert(trialTransformations[i].extent(2) == quadPointCount);
    }
    assert(weights.size() == quadPointCount);

    result.set_size(resultComponentCount, trialDofCount, evalPointCount);
    std::fill(result.begin(), result.end(), 0.);

    std::vector<ResultType> value(resultComponentCount);

    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof)
        for (size_t evalPoint = 0; evalPoint < evalPointCount; ++evalPoint)
            for (size_t quadPoint = 0; quadPoint < quadPointCount; ++quadPoint) {
                m_functor.evaluate(
                            trialGeomData.const_slice(quadPoint),
                            kernels.const_slice(evalPoint, quadPoint),
                            trialTransformations.const_slice(trialDof, quadPoint),
                            value);
                for (int dim = 0; dim < resultComponentCount; ++dim)
                    result(dim, trialDof, evalPoint) += value[dim] *
                            trialGeomData.integrationElements(quadPoint) *
                            weights[quadPoint];
            }
}

} // namespace Fiber

#endif
