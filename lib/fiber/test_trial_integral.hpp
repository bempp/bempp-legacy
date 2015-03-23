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

#ifndef fiber_test_trial_integral_hpp
#define fiber_test_trial_integral_hpp

#include "../common/common.hpp"
#include "types.hpp"

#include "scalar_traits.hpp"

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename T> class CollectionOf3dArrays;
template <typename CoordinateType> class GeometricalData;
/** \endcond */

template <typename BasisFunctionType_, typename ResultType_>
class TestTrialIntegral {
public:
  typedef BasisFunctionType_ BasisFunctionType;
  typedef ResultType_ ResultType;
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  virtual ~TestTrialIntegral() {}

  virtual void addGeometricalDependencies(size_t &geomDeps) const = 0;

  virtual void
  evaluate(const GeometricalData<CoordinateType> &geomData,
           const CollectionOf3dArrays<BasisFunctionType> &testValues,
           const CollectionOf3dArrays<BasisFunctionType> &trialValues,
           const std::vector<CoordinateType> &weights,
           Matrix<ResultType> &result) const = 0;
};

} // namespace Fiber

#endif
