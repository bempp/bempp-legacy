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

#include "../common/common.hpp"

#include "scalar_traits.hpp"

#include "../common/armadillo_fwd.hpp"
#include <vector>

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename T> class CollectionOf3dArrays;
template <typename T> class CollectionOf4dArrays;
template <typename CoordinateType> class GeometricalData;
/** \endcond */

/** \ingroup weak_form_elements
 *  \brief An integral representing the weak form of an integral operator.
 *
 *  This class represents the integral
 *
 *  \f[ \int_\Gamma \int_\Sigma I(x, y)\, d\Gamma(x)\, d\Sigma(y), \f]
 *
 *  defined on a pair of test and trial elements \f$(\Gamma, \Sigma)\f$,
 *  where the integrand \f$f(x, y)\f$ may depend on any number of kernels, test
 *  and trial function transformations, and on any geometrical data
 *  related to the points \f$x\f$ and \f$y\f$ (e.g. their global coordinates or
 *  the unit vector normal to \f$\Gamma\f$ and \f$\Sigma\f$ at these points).
 *
 *  \tparam BasisFunctionType_
 *    Type of the values of the (components of the) basis functions occurring
 *    in the integral (possibly in a transformed form).
 *  \tparam KernelType_
 *    Type of the values of the (components of the) kernel functions occurring
 *    in the integral.
 *  \tparam ResultType_
 *    Type used to represent the value of the integral.
 *
 *  All three template parameters can take the following values: \c float, \c
 *  double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>. All
 *  types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. If either \p
 *  BasisFunctionType_ or \p KernelType_ is a complex type, then \p ResultType_
 *  must be set to the same type. */
template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class TestKernelTrialIntegral {
public:
  typedef BasisFunctionType_ BasisFunctionType;
  typedef KernelType_ KernelType;
  typedef ResultType_ ResultType;
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  /** \brief Destructor. */
  virtual ~TestKernelTrialIntegral() {}

  /** \brief Retrieve types of geometrical data on which the integrand of this
   *  integral depends explicitly.
   *
   *  An implementation of this function for a particular integral
   *  should modify the \p testGeomDeps and \p trialGeomDeps bitfields by
   *  adding to them, using the bitwise OR operation, an appropriate
   *  combination of the flags defined in the enum GeometricalDataType.
   *
   *  For example, an integral whose integrand depends explicitly on the global
   *  coordinates of test and trial points and on the orientation of the
   *  vector normal to the trial element at trial points should modify the
   *  arguments as follows:
      \code
      testGeomDeps |= GLOBALS;
      trialGeomDeps |= GLOBALS | NORMALS;
      \endcode

      \note It is only necessary to modify \p testGeomDeps and \p
      trialGeomDeps if the geometric quantities occur in the integral
      <em>outside any kernels or shape function transformations</em>. For
      example, it is not necessary to add \c NORMALS to \p trialGeomDeps just
      because a weak form contains the double-layer-potential boundary
      operator, which requires the normal to the trial element. It is not
      necessary, either, to ever add the flag INTEGRATION_ELEMENTS as
      integration elements (see Bempp::Geometry::getIntegrationElements() for
      their definition) are automatically included in the list of geometrical
      data required by integrals.
   */
  virtual void addGeometricalDependencies(size_t &testGeomDeps,
                                          size_t &trialGeomDeps) const = 0;

  /** \brief Evaluate the integral using a tensor-product quadrature rule.
   *
   *  This function should evaluate the integral using a quadrature rule of the
   *  form
   *  \f[ \int_\Gamma \int_\Sigma I(x, y)\, d\Gamma(x)\, d\Sigma(y) =
   *      \int_{\hat\Gamma} \int_{\hat\Sigma} I(x(\hat x), y(\hat y)) \,
   *      \mu(\hat x) \, \nu(\hat y) \,
   *      d\hat\Gamma(\hat x)\, d\hat\Sigma(\hat y) \approx
   *      \sum_{p=1}^P \sum_{q=1}^Q w_p w_q \, I(x(\hat x_p), y(\hat y_q)) \,
   *      \mu(\hat x_p) \, \nu(\hat y_q), \f]
   *  where \f$\Gamma\f$ and \f$\Sigma\f$ are the test and trial elements,
   *  \f$x\f$ and \f$y\f$ the physical coordinates on these elements (global
   *  coordinates), \f$\hat\Gamma\f$ and \f$\hat\Sigma\f$ the reference test
   *  and trial elements, \f$\hat x\f$ and \f$\hat y\f$ the coordinates on
   *  the reference elements (local coordinates), \f$\hat x_p\f$ and \f$\hat
   *  y_q\f$ the local coordinates of quadrature points on the test and trial
   *  elements, \f$\hat w_p\f$ and \f$\hat w_q\f$ the corresponding
   *  quadrature weights, and \f$\mu(\hat x_p)\f$ and \f$\nu(\hat y_q)\f$ the
   *  "integration elements" at the quadrature points, defined as
   *  \f$\sqrt{\lvert\det J^T J\rvert}\f$, where \f$J\f$ is the Jacobian matrix
   *  of the local-to-global coordinate mapping.
   *
   *  \param[in] testGeomData
   *    Geometrical data related to the quadrature points on the test element.
   *    The set of available geometrical data always includes integration
   *    elements.
   *  \param[in] trialGeomData
   *    Geometrical data related to the quadrature points on the trial element.
   *    The set of available geometrical data always includes integration
   *    elements.
   *  \param[in] testTransformations
   *    Collection of 3D arrays containing the values of test function
   *    transformations at quadrature points. The number
   *    <tt>testTransformations[i](j, k, p)</tt> is the <em>j</em>th
   *    component of the vector being the value of the <em>i</em>th
   *    transformation of the <em>k</em>th test function at the
   *    <em>p</em>th test quadrature point.
   *  \param[in] trialTransformations
   *    Collection of 3D arrays containing the values of trial function
   *    transformations at quadrature points. The number
   *    <tt>trialTransformations[i](j, k, p)</tt> is the <em>j</em>th
   *    component of the vector being the value of the <em>i</em>th
   *    transformation of the <em>k</em>th trial function at the
   *    <em>p</em>th trial quadrature point.
   *  \param[in] kernels
   *    Collection of 4D arrays containing the values of kernels. The number
   *    <tt>kernels[i][(j, k, p, q)</tt> is the (<em>j</em>, <em>k</em>)th
   *    entry in the tensor being the value of the <em>i</em>th kernel at the
   *    <em>p</em>th test point and <em>q</em>th trial point.
   *  \param[in] testQuadWeights
   *    Vector of the quadrature weights corresponding to the quadrature
   *    points on the test elements.
   *  \param[in] trialQuadWeights
   *    Vector of the quadrature weights corresponding to the quadrature
   *    points on the trial elements.
   *  \param[out] result
   *    Two-dimensional array whose (<em>i</em>, <em>j</em>)th element should
   *    contain, on output, the value of the integral involving the
   *    <em>i</em>th test function and <em>j</em>th trial function.
   */
  virtual void evaluateWithTensorQuadratureRule(
      const GeometricalData<CoordinateType> &testGeomData,
      const GeometricalData<CoordinateType> &trialGeomData,
      const CollectionOf3dArrays<BasisFunctionType> &testTransformations,
      const CollectionOf3dArrays<BasisFunctionType> &trialTransformations,
      const CollectionOf4dArrays<KernelType> &kernels,
      const std::vector<CoordinateType> &testQuadWeights,
      const std::vector<CoordinateType> &trialQuadWeights,
      arma::Mat<ResultType> &result) const = 0;

  /** \brief Evaluate the integral using a non-tensor-product quadrature rule.
   *
   *  This function should evaluate the integral using a quadrature rule of the
   *  form
   *  \f[ \int_\Gamma \int_\Sigma I(x, y)\, d\Gamma(x)\, d\Sigma(y) =
   *      \int_{\hat\Gamma} \int_{\hat\Sigma} I(x(\hat x), y(\hat y)) \,
   *      \mu(\hat x) \, \nu(\hat y) \,
   *      d\hat\Gamma(\hat x)\, d\hat\Sigma(\hat y) \approx
   *      \sum_{p=1}^P w_p \, I(x(\hat x_p), y(\hat y_p)) \,
   *      \mu(\hat x_p) \, \nu(\hat y_p), \f]
   *  where \f$\Gamma\f$ and \f$\Sigma\f$ are the test and trial elements,
   *  \f$x\f$ and \f$y\f$ the physical coordinates on these elements (global
   *  coordinates), \f$\hat\Gamma\f$ and \f$\hat\Sigma\f$ the reference test
   *  and trial elements, \f$\hat x\f$ and \f$\hat y\f$ the coordinates on
   *  the reference elements (local coordinates), \f$\hat x_p\f$ and \f$\hat
   *  y_p\f$ the local coordinates of quadrature points on the test and trial
   *  elements, \f$\hat w_p\f$ the corresponding
   *  quadrature weights, and \f$\mu(\hat x_p)\f$ and \f$\nu(\hat y_p)\f$ the
   *  "integration elements" at the quadrature points, defined as
   *  \f$\sqrt{\lvert\det J^T J\rvert}\f$, where \f$J\f$ is the Jacobian matrix
   *  of the local-to-global coordinate mapping.
   *
   *  \param[in] testGeomData
   *    Geometrical data related to the quadrature points on the test element.
   *    The set of available geometrical data always includes integration
   *    elements.
   *  \param[in] trialGeomData
   *    Geometrical data related to the quadrature points on the trial element.
   *    The set of available geometrical data always includes integration
   *    elements.
   *  \param[in] testTransformations
   *    Collection of 3D arrays containing the values of test function
   *    transformations at quadrature points. The number
   *    <tt>testTransformations[i](j, k, p)</tt> is the <em>j</em>th
   *    component of the vector being the value of the <em>i</em>th
   *    transformation of the <em>k</em>th test function at the
   *    <em>p</em>th test quadrature point.
   *  \param[in] trialTransformations
   *    Collection of 3D arrays containing the values of trial function
   *    transformations at quadrature points. The number
   *    <tt>trialTransformations[i](j, k, p)</tt> is the <em>j</em>th
   *    component of the vector being the value of the <em>i</em>th
   *    transformation of the <em>k</em>th trial function at the
   *    <em>p</em>th trial quadrature point.
   *  \param[in] kernels
   *    Collection of 3D arrays containing the values of kernels. The number
   *    <tt>kernels[i][(j, k, p)</tt> is the (<em>j</em>, <em>k</em>)th
   *    entry in the tensor being the value of the <em>i</em>th kernel at the
   *    <em>p</em>th test and trial point.
   *  \param[in] testQuadWeights
   *    Vector of the quadrature weights corresponding to the quadrature
   *    points.
   *  \param[out] result
   *    Two-dimensional array whose (<em>i</em>, <em>j</em>)th element should
   *    contain, on output, the value of the integral involving the
   *    <em>i</em>th test function and <em>j</em>th trial function.
   */
  virtual void evaluateWithNontensorQuadratureRule(
      const GeometricalData<CoordinateType> &testGeomData,
      const GeometricalData<CoordinateType> &trialGeomData,
      const CollectionOf3dArrays<BasisFunctionType> &testTransformations,
      const CollectionOf3dArrays<BasisFunctionType> &trialTransformations,
      const CollectionOf3dArrays<KernelType> &kernels,
      const std::vector<CoordinateType> &quadWeights,
      arma::Mat<ResultType> &result) const = 0;
};

} // namespace Fiber

#endif
