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

#ifndef fiber_default_kernel_trial_integral_hpp
#define fiber_default_kernel_trial_integral_hpp

#include "kernel_trial_integral.hpp"

namespace Fiber {

template <typename IntegrandFunctor>
class DefaultKernelTrialIntegral
    : public KernelTrialIntegral<typename IntegrandFunctor::BasisFunctionType,
                                 typename IntegrandFunctor::KernelType,
                                 typename IntegrandFunctor::ResultType> {
  typedef KernelTrialIntegral<typename IntegrandFunctor::BasisFunctionType,
                              typename IntegrandFunctor::KernelType,
                              typename IntegrandFunctor::ResultType> Base;

public:
  typedef typename Base::CoordinateType CoordinateType;
  typedef typename Base::BasisFunctionType BasisFunctionType;
  typedef typename Base::KernelType KernelType;
  typedef typename Base::ResultType ResultType;

  explicit DefaultKernelTrialIntegral(const IntegrandFunctor &functor)
      : m_functor(functor) {}

  virtual int resultDimension() const;

  virtual void addGeometricalDependencies(size_t &trialGeomDeps) const;

  virtual void
  evaluate(const GeometricalData<CoordinateType> &trialGeomData,
           const CollectionOf4dArrays<KernelType> &kernels,
           const CollectionOf2dArrays<ResultType> &weightedTrialTransformations,
           const std::vector<CoordinateType> &weights,
           _2dArray<ResultType> &result) const;

  virtual void evaluateWithPureWeights(
      const GeometricalData<CoordinateType> &trialGeomData,
      const CollectionOf4dArrays<KernelType> &kernels,
      const CollectionOf3dArrays<BasisFunctionType> &trialTransformations,
      const std::vector<CoordinateType> &weights,
      _3dArray<ResultType> &result) const;

private:
  IntegrandFunctor m_functor;
};

} // namespace Fiber

#include "default_kernel_trial_integral_imp.hpp"

#endif
