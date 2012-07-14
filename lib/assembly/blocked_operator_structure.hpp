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

template <typename BasisFunctionType, typename ResultType>
class BlockedOperatorStructure
{
public:
    void setBlock(size_t row, size_t column,
                  const BoundaryOperator<BasisFunctionType, ResultType>& op);
    void resetBlock(size_t row, size_t column);
    void reset();

    BoundaryOperator<BasisFunctionType, ResultType> block(
            size_t row, size_t column) const;
    bool isEmpty(size_t row, size_t column) const;

    size_t rowCount() const;
    size_t columnCount() const;

private:
    typedef std::pair<size_t, size_t> Key; // row, column
    typedef BoundaryOperator<BasisFunctionType, ResultType> BoundaryOp;
    std::map<Key, BoundaryOp> m_blocks;
};

} // namespace Bempp

#endif
