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

#ifndef fiber_test_scalar_kernel_trial_integrand_functor_hpp
#define fiber_test_scalar_kernel_trial_integrand_functor_hpp

#include "../common/common.hpp"

#include <cassert>
#include "collection_of_3d_arrays.hpp"
#include "geometrical_data.hpp"
#include "conjugate.hpp"
#include "scalar_traits.hpp"

namespace Fiber {

/** \ingroup functors
 *  \brief Functor evaluating the integrand of most standard boundary integral
 *operators.
 *
 *  ....... TODO .......
 */

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class TestScalarKernelTrialIntegrandFunctor {
public:
  typedef BasisFunctionType_ BasisFunctionType;
  typedef KernelType_ KernelType;
  typedef ResultType_ ResultType;
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  void addGeometricalDependencies(size_t &testGeomDeps,
                                  size_t &trialGeomDeps) const {
    // do nothing
  }

  template <template <typename T> class CollectionOf2dSlicesOfConstNdArrays>
  ResultType evaluate(
      const ConstGeometricalDataSlice<CoordinateType> & /* testGeomData */,
      const ConstGeometricalDataSlice<CoordinateType> & /* trialGeomData */,
      const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType> &testValues,
      const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType> &trialValues,
      const CollectionOf2dSlicesOfConstNdArrays<KernelType> &kernelValues)
      const {
    // Assert that there is at least one scalar-valued kernel
    assert(kernelValues.size() >= 1);
    assert(kernelValues[0].extent(0) == 1);
    assert(kernelValues[0].extent(1) == 1);

    // Assert that there is at least one test and trial transformation
    // and that their dimensions agree
    assert(testValues.size() >= 1);
    assert(trialValues.size() >= 1);
    const size_t transCount = testValues.size();
    assert(trialValues.size() == transCount);
    if (kernelValues.size() > 1)
      assert(kernelValues.size() == transCount);

    ResultType result = 0.;
    for (size_t transIndex = 0; transIndex < transCount; ++transIndex) {
      const size_t transDim = testValues[transIndex].extent(0);
      assert(trialValues[transIndex].extent(0) == transDim);
      BasisFunctionType dotProduct = 0.;
      for (size_t dim = 0; dim < transDim; ++dim)
        dotProduct += conjugate(testValues[transIndex](dim)) *
                      trialValues[transIndex](dim);
      if (transCount == 1)
        result += dotProduct * kernelValues[0](0, 0);
      else
        result += dotProduct * kernelValues[transIndex](0, 0);
    }
    return result;
  }
};

} // namespace Fiber

#endif
