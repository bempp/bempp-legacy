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

#include "blocked_operator_structure.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
void BlockedOperatorStructure<BasisFunctionType, ResultType>::setBlock(
    size_t row, size_t column,
    const BoundaryOperator<BasisFunctionType, ResultType> &op) {
  if (!op.isInitialized())
    resetBlock(row, column);
  else
    m_blocks[Key(row, column)] = op;
}

template <typename BasisFunctionType, typename ResultType>
void BlockedOperatorStructure<BasisFunctionType, ResultType>::resetBlock(
    size_t row, size_t column) {
  m_blocks.erase(Key(row, column));
}

template <typename BasisFunctionType, typename ResultType>
void BlockedOperatorStructure<BasisFunctionType, ResultType>::reset() {
  m_blocks.clear();
}

template <typename BasisFunctionType, typename ResultType>
BoundaryOperator<BasisFunctionType, ResultType>
BlockedOperatorStructure<BasisFunctionType, ResultType>::block(
    size_t row, size_t column) const {
  typedef typename std::map<Key, BoundaryOp>::const_iterator Iterator;
  Iterator it = m_blocks.find(Key(row, column));
  if (it == m_blocks.end())
    return BoundaryOperator<BasisFunctionType, ResultType>();
  else
    // it->second is the boundary operator
    return it->second;
}

template <typename BasisFunctionType, typename ResultType>
bool BlockedOperatorStructure<BasisFunctionType, ResultType>::isEmpty(
    size_t row, size_t column) const {
  return (m_blocks.find(Key(row, column)) == m_blocks.end());
}

template <typename BasisFunctionType, typename ResultType>
size_t
BlockedOperatorStructure<BasisFunctionType, ResultType>::rowCount() const {
  typedef typename std::map<Key, BoundaryOp>::const_iterator Iterator;
  size_t count = 0;
  for (Iterator it = m_blocks.begin(); it != m_blocks.end(); ++it)
    // it->first is the key, it->second is the boundary operator
    if (it->second.isInitialized())
      count = std::max(count, it->first.first +
                                  1 /* count is one more than max index */);
  return count;
}

template <typename BasisFunctionType, typename ResultType>
size_t
BlockedOperatorStructure<BasisFunctionType, ResultType>::columnCount() const {
  typedef typename std::map<Key, BoundaryOp>::const_iterator Iterator;
  size_t count = 0;
  for (Iterator it = m_blocks.begin(); it != m_blocks.end(); ++it)
    // it->first is the key, it->second is the boundary operator
    if (it->second.isInitialized())
      count = std::max(count, it->first.second +
                                  1 /* count is one more than max index */);
  return count;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(BlockedOperatorStructure);

} // namespace Bempp
