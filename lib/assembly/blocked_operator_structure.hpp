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

#ifndef bempp_blocked_operator_structure_hpp
#define bempp_blocked_operator_structure_hpp

#include "../common/common.hpp"

#include "boundary_operator.hpp"

#include <map>

namespace Bempp
{

/** \ingroup boundary_operators
 *  \brief Helper class used in construction of blocked boundary operators.
 *
 *  This class represents a matrix of boundary operators. To construct a
 *  blocked boundary operator, store individual boundary operators in
 *  appropriate rows and columns of a BlockedOperatorStructure object and pass
 *  this object to the constructor of the BlockedBoundaryOperator class for
 *  validation.
 */
template <typename BasisFunctionType, typename ResultType>
class BlockedOperatorStructure
{
public:
    /** \brief Store a boundary operator in a given block.
     *
     *  This function stores operator \p op in the block from row \p row and
     *  column \p column, replacing any operator that might have been stored
     *  there before.
     *
     *  \note Row and column numbers start from 0. */
    void setBlock(size_t row, size_t column,
                  const BoundaryOperator<BasisFunctionType, ResultType>& op);

    /** \brief Reset a block.
     *
     *  This function restores the original (uninitialized) state of the block
     *  from row \p row and column \p column. */
    void resetBlock(size_t row, size_t column);

    /** \brief Reset all blocks.
     *
     *  This function restores the original (uninitialized) state of all
     *  blocks. */
    void reset();

    /** \brief Return the boundary operator stored in a given block.
     *
     *  This function returns the boundary operator stored the block from row
     *  \p row and column \p column, or a null (uninitialized) operator if no
     *  operator has previously been assigned to this block. */
    BoundaryOperator<BasisFunctionType, ResultType> block(
            size_t row, size_t column) const;

    /** \brief Return whether a block is empty.
     *
     *  This function returns true if no operator is stored in the block from
     *  row \p row and column \p column, false otherwise. */
    bool isEmpty(size_t row, size_t column) const;

    /** \brief Return the number of rows in the matrix.
     *
     *  This function returns the smallest number \f$m \geq 0\f$ such that
     *  all the blocks in rows with indices \f$m, m+1, \dots\f$ are empty. */
    size_t rowCount() const;

    /** \brief Return the number of columns in the matrix.
     *
     *  This function returns the smallest number \f$n \geq 0\f$ such that
     *  all the blocks in columns with indices \f$n, n+1, \dots\f$ are empty. */
    size_t columnCount() const;

private:
    /** \cond PRIVATE */
    typedef std::pair<size_t, size_t> Key; // row, column
    typedef BoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;
    std::map<Key, BoundaryOp> m_blocks;
    /** \endcond */
};

} // namespace Bempp

#endif
