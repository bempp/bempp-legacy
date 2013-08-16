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

#ifndef fiber_default_test_kernel_trial_integral_imp_hpp
#define fiber_default_test_kernel_trial_integral_imp_hpp

#include "default_test_kernel_trial_integral.hpp"

#include "geometrical_data.hpp"

#include <cassert>

namespace Fiber
{

template <typename IntegrandFunctor>
void DefaultTestKernelTrialIntegral<IntegrandFunctor>::
addGeometricalDependencies(size_t& testGeomDeps, size_t& trialGeomDeps) const
{
    testGeomDeps |= INTEGRATION_ELEMENTS;
    trialGeomDeps |= INTEGRATION_ELEMENTS;

    m_functor.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
}

template <typename IntegrandFunctor>
void DefaultTestKernelTrialIntegral<IntegrandFunctor>::
evaluateWithTensorQuadratureRule(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        const CollectionOf3dArrays<BasisFunctionType>& testValues,
        const CollectionOf3dArrays<BasisFunctionType>& trialValues,
        const CollectionOf4dArrays<KernelType>& kernelValues,
        const std::vector<CoordinateType>& testQuadWeights,
        const std::vector<CoordinateType>& trialQuadWeights,
        arma::Mat<ResultType>& result) const
{
    // Evaluate constants

    const size_t testDofCount = testValues[0].extent(1);
    const size_t trialDofCount = trialValues[0].extent(1);

    const size_t testPointCount = testQuadWeights.size();
    const size_t trialPointCount = trialQuadWeights.size();

    // Assert that array dimensions are correct

    for (size_t i = 0; i < kernelValues.size(); ++i) {
        assert(kernelValues[i].extent(2) == testPointCount);
        assert(kernelValues[i].extent(3) == trialPointCount);
    }
    for (size_t i = 0; i < testValues.size(); ++i)
        assert(testValues[i].extent(2) == testPointCount);
    for (size_t i = 0; i < trialValues.size(); ++i)
        assert(trialValues[i].extent(2) == trialPointCount);

    // Integrate

    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof)
        for (size_t testDof = 0; testDof < testDofCount; ++testDof)
        {
            ResultType sum = 0.;
            for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
                const CoordinateType trialWeight =
                        trialGeomData.integrationElements(trialPoint) *
                        trialQuadWeights[trialPoint];
                ResultType partialSum = 0.;
                for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
                    const CoordinateType testWeight =
                            testGeomData.integrationElements(testPoint) *
                            testQuadWeights[testPoint];
                    partialSum += m_functor.evaluate(
                                testGeomData.const_slice(testPoint),
                                trialGeomData.const_slice(trialPoint),
                                testValues.const_slice(testDof, testPoint),
                                trialValues.const_slice(trialDof, trialPoint),
                                kernelValues.const_slice(testPoint, trialPoint)) *
                            testWeight;
                }
                sum += partialSum * trialWeight;
            }
            result(testDof, trialDof) = sum;
        }
}

template <typename IntegrandFunctor>
void DefaultTestKernelTrialIntegral<IntegrandFunctor>::
evaluateWithNontensorQuadratureRule(
        const GeometricalData<CoordinateType>& testGeomData,
        const GeometricalData<CoordinateType>& trialGeomData,
        CollectionOf3dArrays<BasisFunctionType>& testValues,
        CollectionOf3dArrays<BasisFunctionType>& trialValues,
        CollectionOf3dArrays<KernelType>& kernelValues,
        const std::vector<CoordinateType>& quadWeights,
        arma::Mat<ResultType>& result) const
{
    // Evaluate constants

    const size_t testDofCount = testValues[0].extent(1);
    const size_t trialDofCount = trialValues[0].extent(1);

    const size_t pointCount = quadWeights.size();

    // Assert that array dimensions are correct

    for (size_t i = 0; i < kernelValues.size(); ++i)
        assert(kernelValues[i].extent(2) == pointCount);
    for (size_t i = 0; i < testValues.size(); ++i)
        assert(testValues[i].extent(2) == pointCount);
    for (size_t i = 0; i < trialValues.size(); ++i)
        assert(trialValues[i].extent(2) == pointCount);

    // Integrate

    for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof)
        for (size_t testDof = 0; testDof < testDofCount; ++testDof)
        {
            ResultType sum = 0.;
            for (size_t point = 0; point < pointCount; ++point)
                sum += m_functor.evaluate(
                            testGeomData.const_slice(point),
                            trialGeomData.const_slice(point),
                            testValues.const_slice(testDof, point),
                            trialValues.const_slice(trialDof, point),
                            kernelValues.const_slice(point)) *
                        (testGeomData.integrationElements(point) *
                         trialGeomData.integrationElements(point) *
                         quadWeights[point]);
            result(testDof, trialDof) = sum;
        }
}

} // namespace Fiber

#endif
