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

#ifndef bempp_solver_hpp
#define bempp_solver_hpp

#include "../common/common.hpp"

#include "solution.hpp"
#include "blocked_solution.hpp"

#include <vector>

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType> class BoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class BlockedBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class GridFunction;

template <typename BasisFunctionType, typename ResultType>
class Solver
{
public:
    enum ConvergenceTestMode {
        TEST_CONVERGENCE_IN_DUAL_TO_RANGE,
        TEST_CONVERGENCE_IN_RANGE
    };

    virtual ~Solver();

    Solution<BasisFunctionType, ResultType> solve(
            const GridFunction<BasisFunctionType, ResultType>& rhs) const {
        return solveImplNonblocked(rhs); 
    }
    BlockedSolution<BasisFunctionType, ResultType> solve(
            const std::vector<GridFunction<BasisFunctionType, ResultType> >&
            rhs) const {
        return solveImplBlocked(rhs); 
    }

protected:
    static void checkConsistency(
        const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
        const GridFunction<BasisFunctionType, ResultType>& rhs,
        ConvergenceTestMode mode);
    static void checkConsistency(
        const BlockedBoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
        const std::vector<GridFunction<BasisFunctionType, ResultType> >& rhs,
        ConvergenceTestMode mode);

    static void constructBlockedGridFunction(
        const arma::Col<ResultType>& solution,
        const BlockedBoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
        std::vector<GridFunction<BasisFunctionType, ResultType> >& solutionFunctions);

private:
    virtual Solution<BasisFunctionType, ResultType> solveImplNonblocked(
        const GridFunction<BasisFunctionType, ResultType>& rhs) const = 0;        
    virtual BlockedSolution<BasisFunctionType, ResultType> solveImplBlocked(
            const std::vector<GridFunction<BasisFunctionType, ResultType> >&
            rhs) const = 0;    
};

} // namespace Bempp

#endif
