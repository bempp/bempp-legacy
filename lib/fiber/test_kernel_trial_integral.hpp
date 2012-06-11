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

#ifndef fiber_test_kernel_trial_integral_hpp
#define fiber_test_kernel_trial_integral_hpp

#include "scalar_traits.hpp"

#include <armadillo>
#include <vector>

namespace Fiber
{

template <typename T> class CollectionOf3dArrays;
template <typename T> class CollectionOf4dArrays;
template <typename CoordinateType> class GeometricalData;

/** \brief Evaluation of double integrals depending on test functions,
 *  kernels and trial functions. */
template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class TestKernelTrialIntegral
{
public:
    typedef BasisFunctionType_ BasisFunctionType;
    typedef KernelType_ KernelType;
    typedef ResultType_ ResultType;
    typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

    virtual ~TestKernelTrialIntegral()
    {}

    virtual void addGeometricalDependencies(
            int& testGeomDeps, int& trialGeomDeps) const = 0;

    virtual void evaluateWithTensorQuadratureRule(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            const CollectionOf3dArrays<BasisFunctionType>& testTransformations,
            const CollectionOf3dArrays<BasisFunctionType>& trialTransformations,
            const CollectionOf4dArrays<KernelType>& kernels,
            const std::vector<CoordinateType>& testQuadWeights,
            const std::vector<CoordinateType>& trialQuadWeights,
            arma::Mat<ResultType>& result) const = 0;

    virtual void evaluateWithNontensorQuadratureRule(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            const CollectionOf3dArrays<BasisFunctionType>& testTransformations,
            const CollectionOf3dArrays<BasisFunctionType>& trialTransformations,
            const CollectionOf3dArrays<KernelType>& kernels,
            const std::vector<CoordinateType>& quadWeights,
            arma::Mat<ResultType>& result) const = 0;
};

} // namespace Fiber

#endif
