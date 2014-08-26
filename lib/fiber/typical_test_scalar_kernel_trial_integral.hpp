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

#ifndef fiber_typical_test_scalar_kernel_trial_integral_hpp
#define fiber_typical_test_scalar_kernel_trial_integral_hpp

#include "test_kernel_trial_integral.hpp"

#include "default_test_kernel_trial_integral.hpp"
#include "test_scalar_kernel_trial_integrand_functor.hpp"

#include <tbb/enumerable_thread_specific.h>

namespace Fiber {

template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class TypicalTestScalarKernelTrialIntegralBase
    : public TestKernelTrialIntegral<BasisFunctionType_, KernelType_,
                                     ResultType_> {
  typedef TestKernelTrialIntegral<BasisFunctionType_, KernelType_, ResultType_>
  Base;

public:
  typedef typename Base::CoordinateType CoordinateType;
  typedef typename Base::BasisFunctionType BasisFunctionType;
  typedef typename Base::KernelType KernelType;
  typedef typename Base::ResultType ResultType;

  virtual void addGeometricalDependencies(size_t &testGeomDeps,
                                          size_t &trialGeomDeps) const;
};

/** \ingroup weak_form_elements
  \brief Implementation of the TestKernelTrialIntegral interface for "typical"
  integrals, taking advantage of BLAS during quadrature.

  This class implements the interface defined by TestKernelTrialIntegral,
  assuming that the integral has the form
  \f[ \int_\Gamma \int_\Sigma \sum_{i=1}^n
      \vec \phi_i(x) \cdot K(x, y) \, \vec \psi_i(y)
      \, d\Gamma(x)\, d\Sigma(y) \f]
  or
  \f[ \int_\Gamma \int_\Sigma \sum_{i=1}^n
      \vec \phi_i(x) \cdot K_i(x, y) \, \vec \psi_i(y)
      \, d\Gamma(x)\, d\Sigma(y) \f]
  where \f$\Gamma\f$ is a test element and \f$\Sigma\f$ a trial element,
  \f$\vec \phi_i\f$ and \f$\vec \psi_i\f$ (\f$i = 1, 2, \cdots, n\f$ wih
  \f$n\f$ an integer) are test and trial function transformations, and \f$K(x,
  y)\f$ or \f$K_i(x, y)\f$ ((\f$i = 1, 2, \cdots, n\f$) are *scalar* kernels.

  The integrals are evaluated numerically; BLAS matrix-matrix multiplication
  routines are used to speed up the process.
 */
template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class TypicalTestScalarKernelTrialIntegral
    : public TypicalTestScalarKernelTrialIntegralBase<
          BasisFunctionType_, KernelType_, ResultType_> {
  // should never be instantiated -- only the specializations (below) should
private:
  TypicalTestScalarKernelTrialIntegral();
};

template <typename BasisFunctionType_, typename ResultType_>
class TypicalTestScalarKernelTrialIntegral<
    BasisFunctionType_, BasisFunctionType_,
    ResultType_> : public TypicalTestScalarKernelTrialIntegralBase<BasisFunctionType_,
                                                                   BasisFunctionType_,
                                                                   ResultType_> {
  typedef TypicalTestScalarKernelTrialIntegralBase<
      BasisFunctionType_, BasisFunctionType_, ResultType_> Base;

public:
  typedef typename Base::CoordinateType CoordinateType;
  typedef typename Base::BasisFunctionType BasisFunctionType;
  typedef typename Base::KernelType KernelType;
  typedef typename Base::ResultType ResultType;

  TypicalTestScalarKernelTrialIntegral() {}

  virtual void evaluateWithTensorQuadratureRule(
      const GeometricalData<CoordinateType> &testGeomData,
      const GeometricalData<CoordinateType> &trialGeomData,
      const CollectionOf3dArrays<BasisFunctionType> &testValues,
      const CollectionOf3dArrays<BasisFunctionType> &trialValues,
      const CollectionOf4dArrays<KernelType> &kernelValues,
      const std::vector<CoordinateType> &testQuadWeights,
      const std::vector<CoordinateType> &trialQuadWeights,
      arma::Mat<ResultType> &result) const;

  virtual void evaluateWithNontensorQuadratureRule(
      const GeometricalData<CoordinateType> &testGeomData,
      const GeometricalData<CoordinateType> &trialGeomData,
      const CollectionOf3dArrays<BasisFunctionType> &testValues,
      const CollectionOf3dArrays<BasisFunctionType> &trialValues,
      const CollectionOf3dArrays<KernelType> &kernelValues,
      const std::vector<CoordinateType> &quadWeights,
      arma::Mat<ResultType> &result) const;
};

template <typename CoordinateType_>
class TypicalTestScalarKernelTrialIntegral<
    CoordinateType_, std::complex<CoordinateType_>,
    std::complex<
        CoordinateType_>> : public TypicalTestScalarKernelTrialIntegralBase<CoordinateType_,
                                                                            std::complex<
                                                                                CoordinateType_>,
                                                                            std::complex<
                                                                                CoordinateType_>> {
  typedef TypicalTestScalarKernelTrialIntegralBase<
      CoordinateType_, std::complex<CoordinateType_>,
      std::complex<CoordinateType_>> Base;

public:
  typedef typename Base::CoordinateType CoordinateType;
  typedef typename Base::BasisFunctionType BasisFunctionType;
  typedef typename Base::KernelType KernelType;
  typedef typename Base::ResultType ResultType;

  TypicalTestScalarKernelTrialIntegral() {}

  virtual void evaluateWithTensorQuadratureRule(
      const GeometricalData<CoordinateType> &testGeomData,
      const GeometricalData<CoordinateType> &trialGeomData,
      const CollectionOf3dArrays<BasisFunctionType> &testValues,
      const CollectionOf3dArrays<BasisFunctionType> &trialValues,
      const CollectionOf4dArrays<KernelType> &kernelValues,
      const std::vector<CoordinateType> &testQuadWeights,
      const std::vector<CoordinateType> &trialQuadWeights,
      arma::Mat<ResultType> &result) const;

  virtual void evaluateWithNontensorQuadratureRule(
      const GeometricalData<CoordinateType> &testGeomData,
      const GeometricalData<CoordinateType> &trialGeomData,
      const CollectionOf3dArrays<BasisFunctionType> &testValues,
      const CollectionOf3dArrays<BasisFunctionType> &trialValues,
      const CollectionOf3dArrays<KernelType> &kernelValues,
      const std::vector<CoordinateType> &quadWeights,
      arma::Mat<ResultType> &result) const;
};

template <typename CoordinateType_>
class TypicalTestScalarKernelTrialIntegral<
    std::complex<CoordinateType_>, CoordinateType_,
    std::complex<
        CoordinateType_>> : public TypicalTestScalarKernelTrialIntegralBase<std::
                                                                                complex<
                                                                                    CoordinateType_>,
                                                                            CoordinateType_,
                                                                            std::complex<
                                                                                CoordinateType_>> {
  typedef TypicalTestScalarKernelTrialIntegralBase<
      std::complex<CoordinateType_>, CoordinateType_,
      std::complex<CoordinateType_>> Base;

public:
  typedef typename Base::CoordinateType CoordinateType;
  typedef typename Base::BasisFunctionType BasisFunctionType;
  typedef typename Base::KernelType KernelType;
  typedef typename Base::ResultType ResultType;

  TypicalTestScalarKernelTrialIntegral();

  // This is the "standard" (non-BLAS-based) implementation
  virtual void evaluateWithTensorQuadratureRule(
      const GeometricalData<CoordinateType> &testGeomData,
      const GeometricalData<CoordinateType> &trialGeomData,
      const CollectionOf3dArrays<BasisFunctionType> &testValues,
      const CollectionOf3dArrays<BasisFunctionType> &trialValues,
      const CollectionOf4dArrays<KernelType> &kernelValues,
      const std::vector<CoordinateType> &testQuadWeights,
      const std::vector<CoordinateType> &trialQuadWeights,
      arma::Mat<ResultType> &result) const;

  virtual void evaluateWithNontensorQuadratureRule(
      const GeometricalData<CoordinateType> &testGeomData,
      const GeometricalData<CoordinateType> &trialGeomData,
      const CollectionOf3dArrays<BasisFunctionType> &testValues,
      const CollectionOf3dArrays<BasisFunctionType> &trialValues,
      const CollectionOf3dArrays<KernelType> &kernelValues,
      const std::vector<CoordinateType> &quadWeights,
      arma::Mat<ResultType> &result) const;

private:
  DefaultTestKernelTrialIntegral<TestScalarKernelTrialIntegrandFunctor<
      BasisFunctionType, KernelType, ResultType>> m_standardIntegral;
};

} // namespace Fiber

#endif
