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

#include "blocked_boundary_operator.hpp"

#include "blocked_operator_structure.hpp"
#include "context.hpp"
#include "discrete_blocked_boundary_operator.hpp"
#include "grid_function.hpp"
#include "../common/boost_make_shared_fwd.hpp"
#include "../common/to_string.hpp"
#include "../common/eigen_support.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../space/space.hpp"

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
BlockedBoundaryOperator<BasisFunctionType, ResultType>::BlockedBoundaryOperator(
    const BlockedOperatorStructure<BasisFunctionType, ResultType> &structure)
    : m_structure(structure) {
  typedef BoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;

  const size_t INVALID = static_cast<size_t>(-1);
  const size_t rowCount = structure.rowCount();
  const size_t columnCount = structure.columnCount();

  m_domains.resize(columnCount);
  m_ranges.resize(rowCount);
  m_dualsToRanges.resize(rowCount);

  for (size_t col = 0; col < columnCount; ++col) {
    size_t firstNonemptyRow = INVALID;
    for (size_t row = 0; row < rowCount; ++row)
      // Is this block non-empty?
      if (shared_ptr<const Space<BasisFunctionType>> domain =
              structure.block(row, col).domain()) {
        if (!m_domains[col]) {
          // Set this->m_domain to this block's domain
          m_domains[col] = domain;
          firstNonemptyRow = row;
        } else
            // Verify that this block's domain matches this->m_domains[col]
            if (!domain->spaceIsCompatible(*m_domains[col]))
          throw std::invalid_argument(
              "BlockedOperatorStructure::BlockedOperatorStructure(): "
              "domains of operators (" +
              toString(firstNonemptyRow) + ", " + toString(col) + ") and (" +
              toString(row) + ", " + toString(col) + ") differ");
      }
    if (firstNonemptyRow == INVALID)
      throw std::invalid_argument(
          "BlockedOperatorStructure::BlockedOperatorStructure(): "
          "column " +
          toString(col) + " is empty");
  }

  for (size_t row = 0; row < rowCount; ++row) {
    size_t firstNonemptyCol = INVALID;
    for (size_t col = 0; col < columnCount; ++col) {
      // Is the block non-empty?
      BoundaryOp block = structure.block(row, col);
      if (block.isInitialized()) {
        if (!m_ranges[row]) {
          // Set this->m_range to this block's range and
          // this->m_dualToRange to this block's dualToRange and
          m_ranges[row] = block.range();
          m_dualsToRanges[row] = block.dualToRange();
          firstNonemptyCol = col;
        } else
            // Verify that this block's range and dualToRange match
            // this->m_ranges[row] and this->m_dualsToRanges[row]
            if (!block.range()->spaceIsCompatible(*m_ranges[row]))
          throw std::invalid_argument(
              "BlockedOperatorStructure::BlockedOperatorStructure(): "
              "ranges of operators (" +
              toString(row) + ", " + toString(firstNonemptyCol) + ") and (" +
              toString(row) + ", " + toString(col) + ") differ");
        if (!block.dualToRange()->spaceIsCompatible(*m_dualsToRanges[row]))
          throw std::invalid_argument(
              "BlockedOperatorStructure::BlockedOperatorStructure(): "
              "duals to ranges of operators (" +
              toString(row) + ", " + toString(firstNonemptyCol) + ") and (" +
              toString(row) + ", " + toString(col) + ") differ");
      }
    }
    if (firstNonemptyCol == INVALID)
      throw std::invalid_argument(
          "BlockedOperatorStructure::BlockedOperatorStructure(): "
          "row " +
          toString(row) + " is empty");
  }

  // Verify that all spaces are non-null
  for (size_t row = 0; row < rowCount; ++row) {
    assert(m_ranges[row]);
    assert(m_dualsToRanges[row]);
  }
  for (size_t col = 0; col < columnCount; ++col)
    assert(m_domains[col]);
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
BlockedBoundaryOperator<BasisFunctionType, ResultType>::block(
    size_t row, size_t column) const {
  return m_structure.block(row, column);
}

template <typename BasisFunctionType, typename ResultType>
bool BlockedBoundaryOperator<BasisFunctionType, ResultType>::isEmpty(
    size_t row, size_t column) const {
  return m_structure.isEmpty(row, column);
}

template <typename BasisFunctionType, typename ResultType>
size_t
BlockedBoundaryOperator<BasisFunctionType, ResultType>::rowCount() const {
  return m_ranges.size();
}

template <typename BasisFunctionType, typename ResultType>
size_t
BlockedBoundaryOperator<BasisFunctionType, ResultType>::columnCount() const {
  return m_domains.size();
}

template <typename BasisFunctionType, typename ResultType>
size_t BlockedBoundaryOperator<
    BasisFunctionType, ResultType>::totalGlobalDofCountInDomains() const {
  size_t sum = 0;
  for (size_t i = 0; i < m_domains.size(); ++i)
    sum += m_domains[i]->globalDofCount();
  return sum;
}

template <typename BasisFunctionType, typename ResultType>
size_t BlockedBoundaryOperator<
    BasisFunctionType, ResultType>::totalGlobalDofCountInRanges() const {
  size_t sum = 0;
  for (size_t i = 0; i < m_ranges.size(); ++i)
    sum += m_ranges[i]->globalDofCount();
  return sum;
}

template <typename BasisFunctionType, typename ResultType>
size_t BlockedBoundaryOperator<
    BasisFunctionType, ResultType>::totalGlobalDofCountInDualsToRanges() const {
  size_t sum = 0;
  for (size_t i = 0; i < m_dualsToRanges.size(); ++i)
    sum += m_dualsToRanges[i]->globalDofCount();
  return sum;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>>
BlockedBoundaryOperator<BasisFunctionType, ResultType>::weakForm() const {
  if (!m_weakForm) {
    m_weakForm = constructWeakForm();
    assert(m_weakForm);
  }
  return m_weakForm;
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType>>
BlockedBoundaryOperator<BasisFunctionType, ResultType>::domain(size_t column)
    const {
  if (column >= m_domains.size())
    throw std::out_of_range(
        "BlockedBoundaryOperator::domain(): invalid column");
  return m_domains[column];
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType>>
BlockedBoundaryOperator<BasisFunctionType, ResultType>::range(size_t row)
    const {
  if (row >= m_ranges.size())
    throw std::out_of_range("BlockedBoundaryOperator::range(): invalid row");
  return m_ranges[row];
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const Space<BasisFunctionType>>
BlockedBoundaryOperator<BasisFunctionType, ResultType>::dualToRange(size_t row)
    const {
  if (row >= m_dualsToRanges.size())
    throw std::out_of_range(
        "BlockedBoundaryOperator::dualToRange(): invalid row");
  return m_dualsToRanges[row];
}

template <typename BasisFunctionType, typename ResultType>
void BlockedBoundaryOperator<BasisFunctionType, ResultType>::apply(
    const TranspositionMode trans,
    const std::vector<GridFunction<BasisFunctionType, ResultType>> &x_in,
    std::vector<GridFunction<BasisFunctionType, ResultType>> &y_inout,
    ResultType alpha, ResultType beta) const {
  const size_t rowCount = this->rowCount();
  const size_t columnCount = this->columnCount();

  // Sanity test
  if (x_in.size() != columnCount)
    throw std::invalid_argument("BlockedBoundaryOperator::apply(): "
                                "x_in has incorrect length");
  if (y_inout.size() != rowCount)
    throw std::invalid_argument("BlockedBoundaryOperator::apply(): "
                                "y_inout has incorrect length");

  for (size_t col = 0; col < columnCount; ++col)
    if (x_in[col].space() != m_domains[col])
      throw std::runtime_error("BlockedBoundaryOperator::apply(): "
                               "grid function x_in[" +
                               toString(col) +
                               "] is incompatible with this operator");

  // TODO: allow uninitialized grid functions
  for (size_t row = 0; row < rowCount; ++row)
    if (y_inout[row].space() != m_ranges[row])
      throw std::runtime_error("BlockedBoundaryOperator::apply(): "
                               "grid function y_inout[" +
                               toString(row) +
                               "] is incompatible with this operator");

  // Build the x vector
  size_t xSize = 0;
  for (size_t col = 0; col < columnCount; ++col)
    xSize += m_domains[col]->globalDofCount();
  Matrix<ResultType> xVals(xSize,1);
  for (size_t col = 0, start = 0; col < columnCount; ++col) {
    const Vector<ResultType> &chunk = x_in[col].coefficients();
    size_t chunkSize = chunk.rows();
    xVals.block(start,0, chunkSize,1 ) = chunk;
    start += chunkSize;
  }

  // Build the y vector
  size_t ySize = 0;
  for (size_t row = 0; row < rowCount; ++row)
    ySize += m_dualsToRanges[row]->globalDofCount();
  Matrix<ResultType> yVals(ySize,1);
  for (size_t row = 0, start = 0; row < rowCount; ++row) {
    Vector<ResultType> chunk =
        y_inout[row].projections(m_dualsToRanges[row]);
    size_t chunkSize = chunk.rows();
    yVals.block(start,0, chunkSize,1 ) = chunk;
    start += chunkSize;
  }

  // Apply operator and assign the result to y_inout's projections

  weakForm()->apply(trans, xVals, yVals, alpha, beta);

  // Assign the result to the grid functions from y_inout
  for (size_t row = 0, start = 0; row < rowCount; ++row) {
    size_t chunkSize = m_dualsToRanges[row]->globalDofCount();
    y_inout[row].setProjections(m_dualsToRanges[row],
                                yVals.col(0).segment(start,chunkSize));
    start += chunkSize;
  }
}

template <typename BasisFunctionType, typename ResultType>
shared_ptr<const DiscreteBoundaryOperator<ResultType>>
BlockedBoundaryOperator<BasisFunctionType, ResultType>::constructWeakForm()
    const {
  typedef DiscreteBoundaryOperator<ResultType> DiscreteOp;
  typedef BoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;

  const size_t rowCount = this->rowCount();
  const size_t columnCount = this->columnCount();
  Fiber::_2dArray<shared_ptr<const DiscreteOp>> blocks(rowCount, columnCount);
  for (size_t col = 0; col < columnCount; ++col)
    for (size_t row = 0; row < rowCount; ++row) {
      BoundaryOp op = m_structure.block(row, col);
      if (op.isInitialized())
        blocks(row, col) = op.weakForm();
    }

  std::vector<size_t> rowCounts(rowCount);
  for (size_t row = 0; row < rowCount; ++row)
    rowCounts[row] = m_dualsToRanges[row]->globalDofCount();
  std::vector<size_t> columnCounts(columnCount);
  for (size_t col = 0; col < columnCount; ++col)
    columnCounts[col] = m_domains[col]->globalDofCount();
  return boost::make_shared<DiscreteBlockedBoundaryOperator<ResultType>>(
      blocks, rowCounts, columnCounts);
}

template <typename BasisFunctionType, typename ResultType>
std::vector<GridFunction<BasisFunctionType, ResultType>> operator*(
    const BlockedBoundaryOperator<BasisFunctionType, ResultType> &op,
    const std::vector<GridFunction<BasisFunctionType, ResultType>> &funs) {
  if (funs.size() != op.columnCount())
    throw std::invalid_argument(
        "operator*(BlockedBoundaryOperator, GridFunction): "
        "'funs' has incorrect length");
  for (size_t i = 0; i < funs.size(); ++i)
    if (!funs[i].isInitialized())
      throw std::invalid_argument(
          "operator*(BlockedBoundaryOperator, GridFunction): "
          "at least one function in 'funs' is uninitialized");

  typedef GridFunction<BasisFunctionType, ResultType> GF;

  std::vector<GF> result(op.rowCount());
  for (size_t i = 0; i < op.rowCount(); ++i) {
    shared_ptr<const Space<BasisFunctionType>> space = op.range(i);
    Vector<ResultType> coefficients(space->globalDofCount());
    coefficients.setZero();
    result[i] = GF(funs[0].context(), space, coefficients);
  }
  op.apply(NO_TRANSPOSE, funs, result, 1., 0.);
  return result;
}

#define INSTANTIATE_NONMEMBER_FUNCTIONS(BASIS, RESULT)                         \
  template std::vector<GridFunction<BASIS, RESULT>> operator*(                 \
      const BlockedBoundaryOperator<BASIS, RESULT> &,                          \
      const std::vector<GridFunction<BASIS, RESULT>> &)
FIBER_ITERATE_OVER_BASIS_AND_RESULT_TYPES(INSTANTIATE_NONMEMBER_FUNCTIONS);
FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(BlockedBoundaryOperator);

} // namespace Bempp
