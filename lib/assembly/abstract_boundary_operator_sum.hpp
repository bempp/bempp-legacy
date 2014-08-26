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

#ifndef bempp_abstract_boundary_operator_sum_hpp
#define bempp_abstract_boundary_operator_sum_hpp

#include "../common/common.hpp"

#include "abstract_boundary_operator_superposition_base.hpp"
#include "boundary_operator.hpp"

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForOperators;
/** \endcond */

} // namespace Fiber

namespace Bempp {

/** \ingroup composite_boundary_operators
 *  \brief Sum of two abstract boundary operators.
 *
 *  This class represents a sum of two boundary operators. */
template <typename BasisFunctionType_, typename ResultType_>
class AbstractBoundaryOperatorSum
    : public AbstractBoundaryOperatorSuperpositionBase<BasisFunctionType_,
                                                       ResultType_> {
  typedef AbstractBoundaryOperatorSuperpositionBase<BasisFunctionType_,
                                                    ResultType_> Base;

public:
  /** \copydoc AbstractBoundaryOperator::BasisFunctionType */
  typedef typename Base::BasisFunctionType BasisFunctionType;
  /** \copydoc AbstractBoundaryOperator::ResultType */
  typedef typename Base::ResultType ResultType;
  /** \copydoc AbstractBoundaryOperator::CoordinateType */
  typedef typename Base::CoordinateType CoordinateType;
  /** \copydoc AbstractBoundaryOperator::QuadratureStrategy */
  typedef typename Base::QuadratureStrategy QuadratureStrategy;
  /** \copydoc AbstractBoundaryOperatorSuperpositionBase::LocalAssembler */
  typedef typename Base::LocalAssembler LocalAssembler;

  /** \brief Constructor.
   *
   *  Construct an operator representing the sum \f$M \equiv L_1 + L_2\f$
   *  of two boundary operators \f$L_1 : X \to Y\f$ and \f$L_2 : X
   *  \to Y\f$.
   *
   *  \param[in] term1 Operator \f$L_1\f$.
   *  \param[in] term2 Operator \f$L_2\f$.
   *  \param[in] symmetry
   *    (Optional) Symmetry of the weak form of the composite operator.
   *    By default is taken as the logical product of the symmetries of the
   *    two operands. Can be set to any combination of the flags defined in
   *    the enumeration type Symmetry.
   *
   *  \note Both terms must be initialized and their domains, ranges and
   *  spaces dual to ranges must be identical, otherwise an exception is
   *  thrown. */
  AbstractBoundaryOperatorSum(
      const BoundaryOperator<BasisFunctionType, ResultType> &term1,
      const BoundaryOperator<BasisFunctionType, ResultType> &term2,
      int symmetry = AUTO_SYMMETRY);

  virtual bool isLocal() const;

  BoundaryOperator<BasisFunctionType_, ResultType_> term1() const;
  BoundaryOperator<BasisFunctionType_, ResultType_> term2() const;

private:
  BoundaryOperator<BasisFunctionType, ResultType> m_term1, m_term2;
};

} // namespace Bempp

#endif
