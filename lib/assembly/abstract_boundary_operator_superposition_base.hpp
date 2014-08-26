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

#ifndef bempp_abstract_boundary_operator_superposition_base_hpp
#define bempp_abstract_boundary_operator_superposition_base_hpp

#include "../common/common.hpp"

#include "abstract_boundary_operator.hpp"
#include "boundary_operator.hpp"

namespace Fiber {

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForIntegralOperators;
/** \endcond */

} // namespace Fiber

namespace Bempp {

/** \ingroup composite_boundary_operators
 *  \brief Base class for abstract boundary operator superpositions.
 *
 *  This class serves as a base class for AbstractBoundaryOperatorSum
 *  and ScaledAbstractBoundaryOperator. It is able to assemble both separated
 *  and joint versions of superpositions of operators. */
template <typename BasisFunctionType_, typename ResultType_>
class AbstractBoundaryOperatorSuperpositionBase
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
   * Fiber::LocalAssemblerForOperators. */
  typedef typename Fiber::LocalAssemblerForIntegralOperators<ResultType>
  LocalAssembler;

  /** \copydoc AbstractBoundaryOperator::AbstractBoundaryOperator */
  AbstractBoundaryOperatorSuperpositionBase(
      const shared_ptr<const Space<BasisFunctionType>> &domain,
      const shared_ptr<const Space<BasisFunctionType>> &range,
      const shared_ptr<const Space<BasisFunctionType>> &dualToRange,
      const std::string &label, int symmetry);

protected:
  virtual shared_ptr<DiscreteBoundaryOperator<ResultType_>>
  assembleWeakFormImpl(const Context<BasisFunctionType, ResultType> &context)
      const;

private:
  shared_ptr<DiscreteBoundaryOperator<ResultType_>>
  assembleJointOperatorWeakFormInDenseMode(
      std::vector<BoundaryOperator<BasisFunctionType, ResultType>> &ops,
      std::vector<ResultType> &opWeights, bool verbose) const;
  shared_ptr<DiscreteBoundaryOperator<ResultType_>>
  assembleJointOperatorWeakFormInAcaMode(
      const Context<BasisFunctionType, ResultType> &context,
      std::vector<BoundaryOperator<BasisFunctionType, ResultType>> &ops,
      std::vector<ResultType> &opWeights) const;
};

} // namespace Bempp

#endif
