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

#ifndef fiber_default_test_kernel_trial_integral_hpp
#define fiber_default_test_kernel_trial_integral_hpp

#include "test_kernel_trial_integral.hpp"

namespace Fiber
{

/** \brief Standard implementation of the TestKernelTrialIntegral interface.

  The <tt>evaluate*()</tt> member functions can be used to evaluate integrals
  of the form

  \f[ \int_\Gamma_i dx \int_\Gamma_j dy \, I[f(x), g(y)] \f],

  where the integrand \f$I(x, y)\f$ depends in an arbitrary way on
  the test function \f$f(x)\f$ and the trial function \f$f(y)\f$.

  This class template takes a single template parameter, \p IntegrandFunctor, which
  should stand for a type providing the following interface:

  \code{.cpp}
class IntegrandFunctor
{
public:
    typedef ... BasisFunctionType;
    typedef ... KernelType;
    typedef ... ResultType;
    typedef ... CoordinateType;

    void addGeometricalDependencies(size_t& testGeomDeps, size_t& trialGeomDeps) const;

    ResultType evaluate(
            const GeometricalDataSlice<CoordinateType>& testGeomData,
            const GeometricalDataSlice<CoordinateType>& trialGeomData,
            const Slice1dCollection<BasisFunctionType>& testTransformations,
            const Slice1dCollection<BasisFunctionType>& trialTransformations,
            const Slice2dCollection<KernelType>& kernels) const;
};
  \endcode

  The addGeometricalDependencies() method should specify any geometrical data
  on which the integral depends explicitly (not through the kernels or basis function
  transformations). For example, if the integrand depends on the vectors normal
  to the surface at test and trial points, the function could have the form

  \code{.cpp}
void addGeometricalDependencies(size_t& testGeomDeps, size_t& trialGeomDeps) const
{
    testGeomDeps |= NORMALS;
    trialGeomDeps |= NORMALS;
}
  \endcode

  The evaluate() method should compute the integrand -- excluding the quadrature
  weights! -- at a single (test point, trial point) pair. It is supplied with the
  following parameters:

  \param[in] testGeomData  Geometric data at a single test point.
  \param[in] trialGeomData Geometric data at a single trial point.
  \param[in] testTransformations
    Values of a collection of transformations of a single test function at the test point.
  \param[in] trialTransformations
    Values of a collection of transformations of a single trial function at the trial point.
  \param[in] kernels
    Values of a collection of kernels at the (test point, trial point) pair.

  \todo Perhaps implement IntegrandFunctor using the curiously recurring template pattern.
 */
template <typename IntegrandFunctor>
class DefaultTestKernelTrialIntegral :
        public TestKernelTrialIntegral<
        typename IntegrandFunctor::BasisFunctionType,
        typename IntegrandFunctor::KernelType,
        typename IntegrandFunctor::ResultType>
{
    typedef TestKernelTrialIntegral<
    typename IntegrandFunctor::BasisFunctionType,
    typename IntegrandFunctor::KernelType,
    typename IntegrandFunctor::ResultType>
    Base;
public:
    typedef typename Base::CoordinateType CoordinateType;
    typedef typename Base::BasisFunctionType BasisFunctionType;
    typedef typename Base::KernelType KernelType;
    typedef typename Base::ResultType ResultType;

    explicit DefaultTestKernelTrialIntegral(const IntegrandFunctor& functor) :
        m_functor(functor)
    {}

    virtual void addGeometricalDependencies(
            size_t& testGeomDeps, size_t& trialGeomDeps) const;

    virtual void evaluateWithTensorQuadratureRule(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            const CollectionOf3dArrays<BasisFunctionType>& testValues,
            const CollectionOf3dArrays<BasisFunctionType>& trialValues,
            const CollectionOf4dArrays<KernelType>& kernelValues,
            const std::vector<CoordinateType>& testQuadWeights,
            const std::vector<CoordinateType>& trialQuadWeights,
            arma::Mat<ResultType>& result) const;

    virtual void evaluateWithNontensorQuadratureRule(
            const GeometricalData<CoordinateType>& testGeomData,
            const GeometricalData<CoordinateType>& trialGeomData,
            const CollectionOf3dArrays<BasisFunctionType>& testValues,
            const CollectionOf3dArrays<BasisFunctionType>& trialValues,
            const CollectionOf3dArrays<KernelType>& kernelValues,
            const std::vector<CoordinateType>& quadWeights,
            arma::Mat<ResultType>& result) const;

private:
    IntegrandFunctor m_functor;
};

} // namespace Fiber

#include "default_test_kernel_trial_integral_imp.hpp"

#endif
