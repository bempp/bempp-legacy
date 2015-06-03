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

#ifndef bempp_elementary_potential_operator_hpp
#define bempp_elementary_potential_operator_hpp

#include "potential_operator.hpp"

#include "../common/shared_ptr.hpp"
#include "../common/eigen_support.hpp"

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename CoordinateType> class CollectionOfBasisTransformations;
template <typename CoordinateType> class CollectionOfShapesetTransformations;
template <typename KernelType> class CollectionOfKernels;
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class KernelTrialIntegral;
template <typename ResultType> class EvaluatorForIntegralOperators;
template <typename ResultType> class LocalAssemblerForPotentialOperators;
/** \endcond */

} // namespace Bempp

namespace Bempp {

/** \cond FORWARD_DECL */
template <typename ValueType> class DiscreteBoundaryOperator;
/** \endcond */

/** \ingroup potential_operators
 *  \brief Elementary potential operator.
 *
 *  This class provides the interface for evaluation of a potential \f$k(x)\f$
 *  defined by the formula
 *
 *  \f[ k(x) = \int_\Gamma F[x, \psi(y)] \, \mathrm{d}\Gamma, \f]
 *
 *  where the integration goes over a surface \f$\Gamma\f$ and the integrand
 *  \f$F\f$ depends on the coordinates of a point \f$x\f$ lying outside
 *  \f$\Gamma\f$ and a surface charge distribution \f$\psi(y)\f$ does not lie on
 *  \f$\Gamma\f$. The function \f$F[x, \psi(y)]\f$ is assumed to be linear in
 *  \f$psi(y)\f$. In the simplest and most common case, \f$F[x, \psi(y)]\f$
 *  is just
 *  \f[
 *    F[x, \psi(y)] = K(x, y) \psi(y),
 *  \f]
 *
 *  where \f$K(x, y)\f$ is a *kernel function*. For more complex operators,
 *  \f$F\f$ might involve some transformations of the charge distribution (e.g.
 *  its surface divergence or curl), the kernel function might be a tensor, or
 *  \f$F\f$ might consist of several terms. The form of \f$F\f$ for a particular
 *  potential operator is determined by the implementation of the virtual
 *  functions integral(), kernels() and trialTransformations().
 *
 *  \tparam BasisFunctionType_
 *    Type of the values of the (components of the) basis functions into
 *    which functions acted upon by the operator are expanded.
 *  \tparam KernelType_
 *    Type of the values of the (components of the) kernel functions occurring
 *    in the integrand of the operator.
 *  \tparam ResultType_
 *    Type of the values of the (components of the) potential.
 *
 *  All three template parameters can take the following values: \c float, \c
 *  double, <tt>std::complex<float></tt> and <tt>std::complex<double></tt>. All
 *  types must have the same precision: for instance, mixing \c float with
 *  <tt>std::complex<double></tt> is not allowed. If either \p
 *  BasisFunctionType_ or \p KernelType_ is a complex type, then \p ResultType_
 *  must be set to the same type. */
template <typename BasisFunctionType_, typename KernelType_,
          typename ResultType_>
class ElementaryPotentialOperator
    : public PotentialOperator<BasisFunctionType_, ResultType_> {
  typedef PotentialOperator<BasisFunctionType_, ResultType_> Base;

public:
  /** \copydoc PotentialOperator::BasisFunctionType */
  typedef typename Base::BasisFunctionType BasisFunctionType;
  /** \brief Type of the values of the (components of the) kernel functions. */
  typedef KernelType_ KernelType;
  /** \copydoc PotentialOperator::ResultType */
  typedef typename Base::ResultType ResultType;
  /** \copydoc PotentialOperator::CoordinateType */
  typedef typename Base::CoordinateType CoordinateType;
  /** \copydoc PotentialOperator::QuadratureStrategy */
  typedef typename Base::QuadratureStrategy QuadratureStrategy;
  /** \brief Type of the appropriate instantiation of
   *  Fiber::EvaluatorForIntegralOperators. */
  typedef Fiber::EvaluatorForIntegralOperators<ResultType> Evaluator;
  /** \brief Type of the appropriate instantiation of
   *  Fiber::LocalAssemblerForPotentialOperators. */
  typedef Fiber::LocalAssemblerForPotentialOperators<ResultType> LocalAssembler;
  /** \brief Type of the appropriate instantiation of
   * Fiber::CollectionOfShapesetTransformations. */
  typedef Fiber::CollectionOfShapesetTransformations<CoordinateType>
      CollectionOfShapesetTransformations;
  /** \brief Type of the appropriate instantiation of
   *Fiber::CollectionOfBasisTransformations.
   *
   *  \deprecated This type is deprecated; use
   *CollectionOfShapesetTransformations
   *  instead. */
  typedef Fiber::CollectionOfShapesetTransformations<CoordinateType>
      CollectionOfBasisTransformations;
  /** \brief Type of the appropriate instantiation of
   *  Fiber::CollectionOfKernels. */
  typedef Fiber::CollectionOfKernels<KernelType> CollectionOfKernels;
  /** \brief Type of the appropriate instantiation of
   *  Fiber::KernelTrialIntegral. */
  typedef Fiber::KernelTrialIntegral<BasisFunctionType, KernelType, ResultType>
      KernelTrialIntegral;

  virtual std::unique_ptr<InterpolatedFunction<ResultType_>>
  evaluateOnGrid(const GridFunction<BasisFunctionType, ResultType> &argument,
                 const Grid &evaluationGrid,
                 const QuadratureStrategy &quadStrategy,
                 const EvaluationOptions &options) const;

  virtual Matrix<ResultType_>
  evaluateAtPoints(const GridFunction<BasisFunctionType, ResultType> &argument,
                   const Matrix<CoordinateType> &evaluationPoints,
                   const QuadratureStrategy &quadStrategy,
                   const EvaluationOptions &options) const;

  virtual AssembledPotentialOperator<BasisFunctionType_, ResultType_>
  assemble(const shared_ptr<const Space<BasisFunctionType>> &space,
           const shared_ptr<const Matrix<CoordinateType>> &evaluationPoints,
           const QuadratureStrategy &quadStrategy,
           const EvaluationOptions &options) const;

  virtual int componentCount() const;

private:
  /** \brief Return the collection of kernel functions occurring in the
   *  integrand of this operator. */
  virtual const CollectionOfKernels &kernels() const = 0;
  /** \brief Return the collection of transformations of the charge
   *  distribution that occur in the weak form of this operator. */
  virtual const CollectionOfBasisTransformations &
  trialTransformations() const = 0;
  /** \brief Return an object representing the integral used to evaluate the
   *  potential of a charge distribution.
   *
   *  Subclasses of #KernelTrialIntegral implement functions that evaluate
   *  the integral using the data provided by a #CollectionOfKernels
   *  representing the kernel functions occurring in the integrand and a
   *  #CollectionOfBasisTransformations representing the charge-distribution
   *  transformations occurring in the integrand. */
  virtual const KernelTrialIntegral &integral() const = 0;

  /** \cond PRIVATE */
  std::unique_ptr<Evaluator>
  makeEvaluator(const GridFunction<BasisFunctionType, ResultType> &argument,
                const QuadratureStrategy &quadStrategy,
                const EvaluationOptions &options) const;

  std::unique_ptr<LocalAssembler>
  makeAssembler(const Space<BasisFunctionType> &space,
                const Matrix<CoordinateType> &evaluationPoints,
                const QuadratureStrategy &quadStrategy,
                const EvaluationOptions &options) const;

  shared_ptr<DiscreteBoundaryOperator<ResultType_>>
  assembleOperator(const Space<BasisFunctionType> &space,
                   const Matrix<CoordinateType> &evaluationPoints,
                   LocalAssembler &assembler,
                   const EvaluationOptions &options) const;

  std::unique_ptr<DiscreteBoundaryOperator<ResultType_>>
  assembleOperatorInDenseMode(const Space<BasisFunctionType> &space,
                              const Matrix<CoordinateType> &evaluationPoints,
                              LocalAssembler &assembler,
                              const EvaluationOptions &options) const;

  std::unique_ptr<DiscreteBoundaryOperator<ResultType_>>
  assembleOperatorInAcaMode(const Space<BasisFunctionType> &space,
                            const Matrix<CoordinateType> &evaluationPoints,
                            LocalAssembler &assembler,
                            const EvaluationOptions &options) const;
  /** \endcond */
};

} // namespace Bempp

#endif
