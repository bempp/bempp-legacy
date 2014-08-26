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

#ifndef bempp_general_elementary_singular_integral_operator_hpp
#define bempp_general_elementary_singular_integral_operator_hpp

#include "elementary_singular_integral_operator.hpp"

namespace Bempp {

/** \ingroup abstract_boundary_operators
 *  \brief Standard implementation of an elementary singular integral operator.
 *
 *  This class provides an implementation of the interface defined by
 *  ElementarySingularIntegralOperator that is sufficient for most purposes. The
 *  constructor takes four functor objects representing the four elements of the
 *  operator's weak form (collection of kernels, collections of test and trial
 *  function transformations, and the weak form integrand). These functors
 *  are used to construct instances of appropriate instantiations of
 *  DefaultCollectionOfKernels, DefaultCollectionOfShapesetTransformations and
 *  DefaultTestKernelTrialIntegral. These objects are stored as private member
 *  variables and are returned by the implementations of the virtual methods
 *  kernels(), testTransformations(), trialTransformations() and integral().
 *
 *  \tparam BasisFunctionType_
 *    Type of the values of the (components of the) basis functions into
 *    which functions acted upon by the operator are expanded.
 *  \tparam KernelType_
 *    Type of the values of the (components of the) kernel functions occurring
 *    in the integrand of the operator.
 *  \tparam ResultType_
 *    Type used to represent elements of the weak form of the operator.
 *
 *  All three template parameters can take the following values: \c float, \c
 *  double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>. All
 *  types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. If either \p
 *  BasisFunctionType_ or \p KernelType_ is a complex type, then \p ResultType_
 *  must be set to the same type. */
template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class GeneralElementarySingularIntegralOperator
    : public ElementarySingularIntegralOperator<BasisFunctionType_, KernelType_,
                                                ResultType_> {
  typedef ElementarySingularIntegralOperator<BasisFunctionType_, KernelType_,
                                             ResultType_> Base;

public:
  /** \brief Type of the values of the basis functions into which functions
   *  acted upon by the operator are expanded. */
  typedef typename Base::BasisFunctionType BasisFunctionType;
  /** \brief Type of the values of kernel functions. */
  typedef typename Base::KernelType KernelType;
  /** \copydoc ElementaryIntegralOperator::ResultType */
  typedef typename Base::ResultType ResultType;
  /** \copydoc ElementaryIntegralOperator::CoordinateType */
  typedef typename Base::CoordinateType CoordinateType;
  /** \copydoc ElementaryIntegralOperator::CollectionOfShapesetTransformations
   */
  typedef typename Base::CollectionOfShapesetTransformations
  CollectionOfShapesetTransformations;
  /** \copydoc ElementaryIntegralOperator::CollectionOfBasisTransformations */
  typedef typename Base::CollectionOfBasisTransformations
  CollectionOfBasisTransformations;
  /** \copydoc ElementaryIntegralOperator::CollectionOfKernels */
  typedef typename Base::CollectionOfKernels CollectionOfKernels;
  /** \copydoc ElementaryIntegralOperator::TestKernelTrialIntegral */
  typedef typename Base::TestKernelTrialIntegral TestKernelTrialIntegral;

  /** \brief Constructor

   *  \param[in] domain
   *    %Function space being the domain of the operator.
   *  \param[in] range
   *    %Function space being the range of the operator.
   *  \param[in] dualToRange
   *    %Function space dual to the the range of the operator.
   *  \param[in] label
   *    Textual label of the operator. If empty, a unique label is generated
   *    automatically.
   *  \param[in] symmetry
   *    Symmetry of the weak form of the operator. Can be any combination of
   *    the flags defined in the enumeration type Symmetry.
   *  \param[in] kernelFunctor
   *    A functor object to be used to evaluate the collection of kernels of
   *    this operator at a single pair of points. The KernelFunctor class
   *    must provide the interface defined in the documentation of
   *    DefaultCollectionOfKernels.
   *  \param[in] testTransformationsFunctor
   *    A functor object to be used to evaluate the collection of test
   *    function transformations at a single point. The
   *    TestTransformationsFunctor class must provide the interface defined in
   *    the documentation of DefaultCollectionOfShapesetTransformations.
   *  \param[in] trialTransformationsFunctor
   *    A functor object to be used to evaluate the collection of trial
   *    function transformations at a single point. The
   *    TrialTransformationsFunctor class must provide the interface defined
   *    in the documentation of DefaultCollectionOfShapesetTransformations.
   *  \param[in] integrandFunctor
   *    A functor object to be used to evaluate the integrand of the weak form
   *    at a single pair of points. The IntegrandFunctor class must provide
   *    the interface defined in the documentation of
   *    DefaultTestKernelTrialIntegral.
   *
   *  None of the shared pointers may be null and the spaces \p range and \p
   *  dualToRange must be defined on the same grid, otherwise an exception is
   *  thrown.
   *
   *  The implementation of this constructor is contained in
   *  general_elementary_singular_integral_operator_imp.hpp. This
   *  header must be included in any file creating a new
   *  GeneralElementarySingularIntegralOperator object.
   */
  template <typename KernelFunctor, typename TestTransformationsFunctor,
            typename TrialTransformationsFunctor, typename IntegrandFunctor>
  GeneralElementarySingularIntegralOperator(
      const shared_ptr<const Space<BasisFunctionType_>> &domain,
      const shared_ptr<const Space<BasisFunctionType_>> &range,
      const shared_ptr<const Space<BasisFunctionType_>> &dualToRange,
      const std::string &label, int symmetry,
      const KernelFunctor &kernelFunctor,
      const TestTransformationsFunctor &testTransformationsFunctor,
      const TrialTransformationsFunctor &trialTransformationsFunctor,
      const IntegrandFunctor &integrandFunctor);

  /** \overload
   *
   *  This constructor takes the same arguments as the preceding one
   *  except for the \p integral argument, which should be a shared pointer
   *  to an instance of (a subclass of) Fiber::TestKernelTrialIntegral.
   */
  template <typename KernelFunctor, typename TestTransformationsFunctor,
            typename TrialTransformationsFunctor>
  GeneralElementarySingularIntegralOperator(
      const shared_ptr<const Space<BasisFunctionType_>> &domain,
      const shared_ptr<const Space<BasisFunctionType_>> &range,
      const shared_ptr<const Space<BasisFunctionType_>> &dualToRange,
      const std::string &label, int symmetry,
      const KernelFunctor &kernelFunctor,
      const TestTransformationsFunctor &testTransformationsFunctor,
      const TrialTransformationsFunctor &trialTransformationsFunctor,
      const shared_ptr<Fiber::TestKernelTrialIntegral<
          BasisFunctionType_, KernelType_, ResultType_>> &integral);

  /** \overload
   *
   *  This constructor takes the same first five arguments as the preceding
   *  ones, but the last four arguments should be shared pointers to
   *  instances of Fiber::CollectionOfKernels,
   *  Fiber::CollectionOfShapesetTransformations and
   *  Fiber::TestKernelTrialIntegral.
   */
  GeneralElementarySingularIntegralOperator(
      const shared_ptr<const Space<BasisFunctionType_>> &domain,
      const shared_ptr<const Space<BasisFunctionType_>> &range,
      const shared_ptr<const Space<BasisFunctionType_>> &dualToRange,
      const std::string &label, int symmetry,
      const shared_ptr<Fiber::CollectionOfKernels<KernelType_>> &kernels,
      const shared_ptr<Fiber::CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<Fiber::CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<Fiber::TestKernelTrialIntegral<
          BasisFunctionType_, KernelType_, ResultType_>> &integral);

  virtual const CollectionOfKernels &kernels() const { return *m_kernels; }
  virtual const CollectionOfShapesetTransformations &
  testTransformations() const {
    return *m_testTransformations;
  }
  virtual const CollectionOfShapesetTransformations &
  trialTransformations() const {
    return *m_trialTransformations;
  }
  virtual const TestKernelTrialIntegral &integral() const {
    return *m_integral;
  }

private:
  /** \cond PRIVATE */
  shared_ptr<CollectionOfKernels> m_kernels;
  shared_ptr<CollectionOfShapesetTransformations> m_testTransformations;
  shared_ptr<CollectionOfShapesetTransformations> m_trialTransformations;
  shared_ptr<TestKernelTrialIntegral> m_integral;
  /** \endcond */
};

} // namespace Bempp

#endif
