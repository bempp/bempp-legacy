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

#include "solver.hpp"

#include "../assembly/grid_function.hpp"
#include "../assembly/boundary_operator.hpp"
#include "../assembly/blocked_boundary_operator.hpp"
#include "../common/to_string.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../space/space.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
Solver<BasisFunctionType, ResultType>::~Solver()
{
}

template <typename BasisFunctionType, typename ResultType>
void Solver<BasisFunctionType, ResultType>::checkConsistency(
        const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
        const GridFunction<BasisFunctionType, ResultType>& rhs,
        ConvergenceTestMode::Mode mode)
{
    if ((mode == ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE &&
         rhs.dualSpace() != boundaryOp.dualToRange()) ||
        (mode == ConvergenceTestMode::TEST_CONVERGENCE_IN_RANGE &&
         rhs.space() != boundaryOp.range()))
        throw std::invalid_argument(
            "Solver::checkConsistency(): spaces do not match");
}

template <typename BasisFunctionType, typename ResultType>
void Solver<BasisFunctionType, ResultType>::checkConsistency(
        const BlockedBoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
        const std::vector<GridFunction<BasisFunctionType, ResultType> >& rhs,
        ConvergenceTestMode::Mode mode)
{
    const size_t columnCount = boundaryOp.columnCount();
    const size_t rowCount = boundaryOp.rowCount();

    if (rhs.size() != rowCount)
        throw std::invalid_argument(
            "Solver::checkConsistency(): incorrect number of grid functions");
    if (mode == ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE)
        ; // don't do anything
    else if (mode == ConvergenceTestMode::TEST_CONVERGENCE_IN_RANGE) {
        for (size_t i = 0; i < rhs.size(); ++i)
            if (rhs[i].space() != boundaryOp.range(i))
                throw std::invalid_argument(
                    "Solver::checkConsistency(): space of grid function #" +
                    toString(i) + 
                    " does not match the range space of the "
                    "corresponding row of the blocked boundary operator");
    }
    else
        throw std::invalid_argument("Invalid convergence testing mode.");
}

template <typename BasisFunctionType, typename ResultType>
std::vector<GridFunction<BasisFunctionType, ResultType> >
Solver<BasisFunctionType, ResultType>::canonicalizeBlockedRhs(
        const BlockedBoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
        const std::vector<GridFunction<BasisFunctionType, ResultType> >& rhs,
        ConvergenceTestMode::Mode mode)
{
    typedef GridFunction<BasisFunctionType, ResultType> GF;

    const size_t functionCount = boundaryOp.rowCount();
    if (rhs.size() != functionCount)
        throw std::invalid_argument("Solver::canonicalizeBlockedRhs(): "
                                    "rhs has incorrect length");

    bool testInDualToRange;
    switch (mode) {
    case ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE:
        testInDualToRange = true;
        break;
    case ConvergenceTestMode::TEST_CONVERGENCE_IN_RANGE:
        testInDualToRange = false;
        break;
    default:
        throw std::invalid_argument("Invalid convergence testing mode.");
    }

    std::vector<GF> result(functionCount);
    for (size_t i = 0; i < functionCount; ++i)
        if (rhs[i].isInitialized()) {
            if (!testInDualToRange && rhs[i].space() != boundaryOp.range(i))
                throw std::invalid_argument(
                        "Solver::canonicalizeBlockedRhs(): space of "
                        "grid function #" + toString(i) +
                        " does not match the range space of the "
                        "corresponding row of the blocked boundary operator");
            result[i] = rhs[i];
        } else { // not initialized
            arma::Col<ResultType> projections(
                            boundaryOp.dualToRange(i)->globalDofCount());
            projections.fill(0.);
            // find an initialized operator in row i
            for (size_t j = 0; j < boundaryOp.columnCount(); ++j) {
                const BoundaryOperator<BasisFunctionType, ResultType>& block =
                        boundaryOp.block(i, j);
                // We rely on the BoundaryOperator's constructor checking
                // that there is a nonempty block in each row
                if (block.isInitialized()) {
                    result[i] = GF(block.context(), block.range(),
                                   block.dualToRange(),
                                   projections);
                    break;
                }
            }
        }
    return result;
}

template <typename BasisFunctionType, typename ResultType>
void Solver<BasisFunctionType, ResultType>::constructBlockedGridFunction(
        const arma::Col<ResultType>& solution,
        const BlockedBoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
        std::vector<GridFunction<BasisFunctionType, ResultType> >& solutionFunctions)
{
    const size_t columnCount = boundaryOp.columnCount();
    const size_t rowCount = boundaryOp.rowCount();

    // Convert chunks of the solution vector into grid functions
    solutionFunctions.resize(columnCount);
    for (size_t i = 0, start = 0; i < solutionFunctions.size(); ++i) {        
        // Find the first non-zero block in column i and retrieve its context
        shared_ptr<const Context<BasisFunctionType, ResultType> > context;
        for (int row = 0; row < rowCount; ++row)
            if (boundaryOp.block(row, i).context()) {
                context = boundaryOp.block(row, i).context();
                break;
            }
        assert(context);

        size_t chunkSize = boundaryOp.domain(i)->globalDofCount();
        solutionFunctions[i] = 
            GridFunction<BasisFunctionType, ResultType>(
                context, 
                boundaryOp.domain(i),
                solution.rows(start, start + chunkSize - 1));
        start += chunkSize;
    }
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Solver);

} // namespace Bempp
