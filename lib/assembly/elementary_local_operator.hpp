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

#ifndef bempp_elementary_local_operator_hpp
#define bempp_elementary_local_operator_hpp

#include "../common/common.hpp"
#include "../common/types.hpp"

#include "abstract_boundary_operator.hpp"

#include "abstract_boundary_operator_id.hpp"
#include <boost/scoped_ptr.hpp>

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename BasisFunctionType> class CollectionOfShapesetTransformations;
template <typename ResultType> class LocalAssemblerForLocalOperators;
template <typename BasisFunctionType, typename ResultType>
class TestTrialIntegral;
/** \endcond */

} // namespace Fiber

namespace Bempp {

/** \cond FORWARD_DECL */
template <typename BasisFunctionType, typename ResultType>
class BoundaryOperator;
/** \endcond */

/** \ingroup identity
 *  \brief Abstract base class of local elementary operators.
 *
 *  See AbstractBoundaryOperator for the documentation of the template
 *  parameters.
 */
template <typename BasisFunctionType_, typename ResultType_>
class ElementaryLocalOperator
    : public AbstractBoundaryOperator<BasisFunctionType_, ResultType_> {
  typedef AbstractBoundaryOperator<BasisFunctionType_, ResultType_> Base;

public:
  /** \copydoc AbstractBoundaryOperator::BasisFunctionType */
  typedef typename Base::BasisFunctionType BasisFunctionType;
  /** \copydoc AbstractBoundaryOperator::ResultType */
  typedef typename Base::ResultType ResultType;
  /** \copydoc AbstractBoundaryOperator::CoordinateType */
  typedef typename Base::CoordinateType CoordinateType;
  /** \copydoc AbstractBoundaryOperator::QuadratureStrategy */
  typedef typename Base::QuadratureStrategy QuadratureStrategy;
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
  /** \brief Type of the appropriate instantiation of Fiber::TestTrialIntegral.
   */
  typedef Fiber::TestTrialIntegral<BasisFunctionType, ResultType>
      TestTrialIntegral;
  /** \brief Type of the appropriate instantiation of
   *  Fiber::LocalAssemblerForLocalOperators. */
  typedef Fiber::LocalAssemblerForLocalOperators<ResultType> LocalAssembler;

  /** \copydoc AbstractBoundaryOperator::AbstractBoundaryOperator */
  ElementaryLocalOperator(
      const shared_ptr<const Space<BasisFunctionType>> &domain,
      const shared_ptr<const Space<BasisFunctionType>> &range,
      const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
      const std::string &label, int symmetry);

  /** \brief Return true. */
  virtual bool isLocal() const;

  /** \brief Create an assembler for the operator. */
  std::unique_ptr<LocalAssembler>
  makeAssembler(const QuadratureStrategy &quadStrategy,
                const AssemblyOptions &options) const;

  /** \brief Overload. */
  std::unique_ptr<LocalAssembler>
  makeAssembler(const ParameterList &parameterList) const;

protected:
  virtual shared_ptr<DiscreteBoundaryOperator<ResultType_>>
  assembleWeakFormImpl(
      const Context<BasisFunctionType, ResultType> &context) const;

private:
  /** \brief Return the collection of test function transformations occurring
   *  in the weak form of this operator. */
  virtual const CollectionOfShapesetTransformations &
  testTransformations() const = 0;

  /** \brief Return the collection of trial function transformations occurring
   *  in the weak form of this operator. */
  virtual const CollectionOfShapesetTransformations &
  trialTransformations() const = 0;

  /** \brief Return an object representing the integral that is the weak form
   *  of this operator.
   *
   *  Subclasses of #TestTrialIntegral implement functions that evaluate
   *  the integral using the data provided by a pair
   *  of #CollectionOfShapesetTransformations objects representing the test and
   *  trial function transformations occurring in the integrand. */
  virtual const TestTrialIntegral &integral() const = 0;

  virtual shared_ptr<DiscreteBoundaryOperator<ResultType_>>
  assembleWeakFormInternalImpl2(
      LocalAssembler &assembler,
      const Context<BasisFunctionType, ResultType> &options) const;

  std::unique_ptr<DiscreteBoundaryOperator<ResultType_>>
  assembleWeakFormInDenseMode(LocalAssembler &assembler,
                              const AssemblyOptions &options) const;

  std::unique_ptr<DiscreteBoundaryOperator<ResultType_>>
  assembleWeakFormInSparseMode(LocalAssembler &assembler,
                               const AssemblyOptions &options) const;

private:
  shared_ptr<const AbstractBoundaryOperatorId> m_id;
};

} // namespace Bempp

#endif
