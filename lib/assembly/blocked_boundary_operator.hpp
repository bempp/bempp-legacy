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

namespace Bempp
{

template <typename BasisFunctionType> class Space;
template <typename ResultType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class BoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class GridFunction;

template <typename BasisFunctionType, typename ResultType>
class BlockedBoundaryOperator
{
public:
    BlockedBoundaryOperator(
            const BlockedOperatorStructure<BasisFunctionType, ResultType>& structure);

    BoundaryOperator<BasisFunctionType, ResultType> block(
            size_t row, size_t column) const;
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

    shared_ptr<const DiscreteBoundaryOperator<ResultType> > weakForm() const;

    shared_ptr<const Space<BasisFunctionType> > domain(size_t column) const;
    shared_ptr<const Space<BasisFunctionType> > range(size_t row) const;
    shared_ptr<const Space<BasisFunctionType> > dualToRange(size_t row) const;

    /** \brief Set <tt>y_inout := alpha * A * x_in + beta * y_inout</tt>, where
     *  \c A is this operator. */
    void apply(const TranspositionMode trans,
               const std::vector<GridFunction<BasisFunctionType, ResultType> >& x_in,
               std::vector<GridFunction<BasisFunctionType, ResultType> >& y_inout,
               ResultType alpha, ResultType beta) const;

private:
    shared_ptr<const DiscreteBoundaryOperator<ResultType> >
    constructWeakForm() const;

private:
    BlockedOperatorStructure<BasisFunctionType, ResultType> m_structure;
    std::vector<shared_ptr<const Space<BasisFunctionType> > > m_domains;
    std::vector<shared_ptr<const Space<BasisFunctionType> > > m_ranges;
    std::vector<shared_ptr<const Space<BasisFunctionType> > > m_dualsToRanges;
    mutable shared_ptr<const DiscreteBoundaryOperator<ResultType> > m_weakForm;
};

} // namespace Bempp

#endif
