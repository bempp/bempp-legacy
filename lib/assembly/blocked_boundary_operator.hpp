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

#ifndef bempp_blocked_boundary_operator_hpp
#define bempp_blocked_boundary_operator_hpp

#include "../common/common.hpp"
#include "../common/shared_ptr.hpp"

#include "blocked_operator_structure.hpp"
#include "transposition_mode.hpp"

#include <vector>

namespace Bempp {

/** \cond FORWARD_DECL */
template <typename BasisFunctionType> class Space;
template <typename ResultType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType, typename ResultType>
class BoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class GridFunction;
/** \endcond */

/** \ingroup boundary_operators
 *  \brief Boundary operator consisting of multiple blocks arranged in a matrix.
 */
template <typename BasisFunctionType, typename ResultType>
class BlockedBoundaryOperator {
public:
  /** \brief Constructor.
   *
   *  \param[in] structure
   *    A BlockedOperatorStructure object determining the operators occupying
   *    individual blocks.
   *
   *  All the boundary operators from a single column of \p structure must
   *  have the same domain, and all the operators from a single row of \p
   *  structure must have the same range and space dual to range. No row and
   *  no column of \p structure must be completely empty (contain only
   *  uninitialized BoundaryOperator objects). If these conditions are not
   *  satisfied, an exception is thrown. */
  BlockedBoundaryOperator(
      const BlockedOperatorStructure<BasisFunctionType, ResultType> &structure);

  /** \brief Return the operator from row \p row and column \p column.
   *
   *  If block (\p row, \p column) is empty, an uninitialized BoundaryOperator
   *  object is returned. */
  BoundaryOperator<BasisFunctionType, ResultType> block(size_t row,
                                                        size_t column) const;
  /** \brief Return whether the block in row \p row and column \p column is
   *  empty. */
  bool isEmpty(size_t row, size_t column) const;

  /** \brief Return number of block rows. */
  size_t rowCount() const;
  /** \brief Return number of block columns. */
  size_t columnCount() const;

  /** \brief Return total number of global degrees of freedom in all domains. */
  size_t totalGlobalDofCountInDomains() const;
  /** \brief Return total number of global degrees of freedom in all ranges. */
  size_t totalGlobalDofCountInRanges() const;
  /** \brief Return total number of global degrees of freedom in all
   *  duals to ranges. */
  size_t totalGlobalDofCountInDualsToRanges() const;

  /** \brief Return the weak form of this boundary operator.
   *
   *  The returned discrete operator represents the matrix
   *  \f[ L =
   *      \begin{bmatrix}
   *        L_{11} & L_{12} & \dots  & L_{1n} \\
   *        L_{21} & L_{22} & \dots  & L_{2n} \\
   *        \vdots & \vdots & \ddots & \vdots \\
   *        L_{m1} & L_{m2} & \dots  & L_{mn}
   *      \end{bmatrix},
   *  \f]
   *  where \f$L_{ij}\f$ is the weak form of the operator from row *i* and
   *  column *j* of this blocked boundary operator. */
  shared_ptr<const DiscreteBoundaryOperator<ResultType>> weakForm() const;

  /** \brief Return the function space being the domain of all the operators
   *  from column \p column of this blocked operator. */
  shared_ptr<const Space<BasisFunctionType>> domain(size_t column) const;
  /** \brief Return the function space being the range of all the operators
   *  from row \p row of this blocked operator. */
  shared_ptr<const Space<BasisFunctionType>> range(size_t row) const;
  /** \brief Return the function space dual to the range of all the operators
   *  from row \p row of this blocked operator. */
  shared_ptr<const Space<BasisFunctionType>> dualToRange(size_t row) const;

  /** \brief Set <tt>y_inout := alpha * A * x_in + beta * y_inout</tt>, where
   *  \c A is this operator. */
  void
  apply(const TranspositionMode trans,
        const std::vector<GridFunction<BasisFunctionType, ResultType>> &x_in,
        std::vector<GridFunction<BasisFunctionType, ResultType>> &y_inout,
        ResultType alpha, ResultType beta) const;

private:
  /** \cond PRIVATE */
  shared_ptr<const DiscreteBoundaryOperator<ResultType>>
  constructWeakForm() const;

private:
  BlockedOperatorStructure<BasisFunctionType, ResultType> m_structure;
  std::vector<shared_ptr<const Space<BasisFunctionType>>> m_domains;
  std::vector<shared_ptr<const Space<BasisFunctionType>>> m_ranges;
  std::vector<shared_ptr<const Space<BasisFunctionType>>> m_dualsToRanges;
  mutable shared_ptr<const DiscreteBoundaryOperator<ResultType>> m_weakForm;
  /** \endcond */
};

/** \relates BlockedBoundaryOperator
 *  \brief Act with a BoundaryOperator on a GridFunction.
 *
 *  This function returns the GridFunction obtained by acting with the operator
 *  \p op on the grid function \p fun. It is equivalent to calling
 *  \code
 *  op.apply(NO_TRANSPOSE, fun, result, 1., 0.);
 *  return result;
 *  \endcode
 *  on GridFunction \p result with space and dual space compatible with
 *  the range and dual to range of \p op. */
template <typename BasisFunctionType, typename ResultType>
std::vector<GridFunction<BasisFunctionType, ResultType>>
operator*(const BlockedBoundaryOperator<BasisFunctionType, ResultType> &op,
          const std::vector<GridFunction<BasisFunctionType, ResultType>> &funs);

} // namespace Bempp

#endif
