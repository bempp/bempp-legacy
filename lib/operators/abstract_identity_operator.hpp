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

#ifndef bempp_abstract_identity_operator_hpp
#define bempp_abstract_identity_operator_hpp

#include "../common/common.hpp"
#include "../common/types.hpp"

#include "../assembly/elementary_local_operator.hpp"

#include "../assembly/abstract_boundary_operator_id.hpp"
#include <boost/scoped_ptr.hpp>

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForOperators;
/** \endcond */

} // namespace Fiber

namespace Bempp {

/** \cond FORWARD_DECL */
template <typename BasisFunctionType, typename ResultType>
class AbstractIdentityOperator;
/** \endcond */

template <typename BasisFunctionType, typename ResultType>
class BEMPP_DEPRECATED AbstractIdentityOperatorId : public AbstractBoundaryOperatorId {
public:
  AbstractIdentityOperatorId(const AbstractIdentityOperator<BasisFunctionType, ResultType> &op);
  virtual size_t hash() const;
  virtual bool isEqual(const AbstractBoundaryOperatorId &other) const;

private:
  const Space<BasisFunctionType> *m_domain;
  const Space<BasisFunctionType> *m_range;
  const Space<BasisFunctionType> *m_dualToRange;
};

/** \ingroup local_operators
 *  \brief Identity operator.
 *
 *  Let \f$X\f$ and \f$Y\f$ be two function spaces defined on the same grid. If
 *  \f$X \supset Y\f$, an instance of AbstractIdentityOperator with domain \f$X\f$
 *  and range \f$Y\f$ represents the orthogonal projection operator from
 *  \f$X\f$ to \f$Y\f$. If \f$X \subset Y\f$, it represents the inclusion
 *  operator from \f$X\f$ to \f$Y\f$. In the special case of \f$X = Y\f$, we
 *  have the standard identity operator. In BEM++ we (ab)use the term "identity
 *  operator" to refer to all these three cases.
 *
 *  See AbstractBoundaryOperator for the documentation of the template
 *  parameters.
 *
 *  Use identityOperator() to create a BoundaryOperator object wrapping
 *  an identity operator.
 */
template <typename BasisFunctionType_, typename ResultType_>
class AbstractIdentityOperator
    : public ElementaryLocalOperator<BasisFunctionType_, ResultType_> {
  typedef ElementaryLocalOperator<BasisFunctionType_, ResultType_> Base;

public:
  /** \copydoc ElementaryLocalOperator::BasisFunctionType */
  typedef typename Base::BasisFunctionType BasisFunctionType;
  /** \copydoc ElementaryLocalOperator::ResultType */
  typedef typename Base::ResultType ResultType;
  /** \copydoc ElementaryLocalOperator::CoordinateType */
  typedef typename Base::CoordinateType CoordinateType;
  /** \copydoc ElementaryLocalOperator::QuadratureStrategy */
  typedef typename Base::QuadratureStrategy QuadratureStrategy;
  /** \copydoc ElementaryLocalOperator::CollectionOfShapesetTransformations */
  typedef typename Base::CollectionOfShapesetTransformations
      CollectionOfShapesetTransformations;
  /** \copydoc ElementaryLocalOperator::CollectionOfBasisTransformations */
  typedef typename Base::CollectionOfBasisTransformations
      CollectionOfBasisTransformations;
  /** \copydoc ElementaryLocalOperator::TestTrialIntegral */
  typedef typename Base::TestTrialIntegral TestTrialIntegral;

  /** \brief Constructor.
   *
   *  \param[in] domain
   *    Function space being the domain of the operator.
   *  \param[in] range
   *    Function space being the range of the operator.
   *  \param[in] dualToRange
   *    Function space dual to the the range of the operator.
   *  \param[in] label
   *    Textual label of the operator. If empty, a unique label is generated
   *    automatically.
   *  \param[in] symmetry
   *    Symmetry of the weak form of the operator. Can be any combination of
   *    the flags defined in the enumeration type Symmetry.
   *    If set to AUTO_SYMMETRY (default), the symmetry is determined
   *    automatically by checking whether its domain and space dual to its
   *    range are equal. If so, the operator is marked as Hermitian,
   *    and if the basis functions are real-valued, also as symmetric.
   *
   *  All the three spaces must be defined on the same grid. */
  AbstractIdentityOperator(
      const shared_ptr<const Space<BasisFunctionType>> &domain,
      const shared_ptr<const Space<BasisFunctionType>> &range,
      const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
      const std::string &label = "", int symmetry = AUTO_SYMMETRY);
  virtual ~AbstractIdentityOperator();

  /** \brief Return the identifier of this operator.
   *
   *  Identity operators are treated as equivalent if they have the same domain,
   *  range and dual to range.
   *
   *  \deprecated This function is deprecated and will be removed in a future
   *  version of BEM++.
   */
  BEMPP_DEPRECATED virtual shared_ptr<const AbstractBoundaryOperatorId>
  id() const;

private:
  virtual const CollectionOfShapesetTransformations &
  testTransformations() const;

  virtual const CollectionOfShapesetTransformations &
  trialTransformations() const;

  virtual const TestTrialIntegral &integral() const;

private:
  shared_ptr<TestTrialIntegral> m_integral;
  shared_ptr<const AbstractBoundaryOperatorId> m_id;
};


} // namespace Bempp

#endif
