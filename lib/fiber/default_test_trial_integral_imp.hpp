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

#ifndef fiber_default_test_trial_integral_imp_hpp
#define fiber_default_test_trial_integral_imp_hpp

#include "default_test_trial_integral.hpp"

#include "geometrical_data.hpp"

#include <algorithm>

namespace Fiber {

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

    void addGeometricalDependencies(int& geomDeps) const;

    ResultType evaluate(
            const ConstGeometricalDataSlice<CoordinateType>& geomData,
            const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>&
            testTransformations,
            const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>&
            trialTransformations) const;
};
*/

template <typename IntegrandFunctor>
void DefaultTestTrialIntegral<IntegrandFunctor>::addGeometricalDependencies(
    size_t &geomDeps) const {
  geomDeps |= INTEGRATION_ELEMENTS;
  m_functor.addGeometricalDependencies(geomDeps);
}

template <typename IntegrandFunctor>
void DefaultTestTrialIntegral<IntegrandFunctor>::evaluate(
    const GeometricalData<CoordinateType> &geomData,
    const CollectionOf3dArrays<BasisFunctionType> &testValues,
    const CollectionOf3dArrays<BasisFunctionType> &trialValues,
    const std::vector<CoordinateType> &weights,
    arma::Mat<ResultType> &result) const {
  const size_t pointCount = weights.size();
  assert(testValues.size() == 1);
  assert(trialValues.size() == 1);

  // We don't care about the number of rows of testValues[0] and trialValues[0]
  // -- it's up to the integrand functor
  const size_t testDofCount = testValues[0].extent(1);
  const size_t trialDofCount = trialValues[0].extent(1);
  assert(testValues[0].extent(2) == pointCount);
  assert(trialValues[0].extent(2) == pointCount);
  for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof)
    for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
      ResultType sum = 0.;
      for (size_t point = 0; point < pointCount; ++point)
        sum += m_functor.evaluate(geomData.const_slice(point),
                                  testValues.const_slice(testDof, point),
                                  trialValues.const_slice(trialDof, point)) *
               geomData.integrationElements(point) * weights[point];
      result(testDof, trialDof) = sum;
    }
}

} // namespace Fiber

#endif
