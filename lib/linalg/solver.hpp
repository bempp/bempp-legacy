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



struct ConvergenceTestMode {
    enum Mode {
        TEST_CONVERGENCE_IN_DUAL_TO_RANGE,
        TEST_CONVERGENCE_IN_RANGE
    };
};

/** \ingroup linalg
  * \brief An abstract interface for various types of solvers
  *
  * This class is an interface to the solution of linear systems in BEM++.
  * Concrete subclasses implement specific linear solvers.
  *
  */

template <typename BasisFunctionType, typename ResultType>
class Solver
{
public:

    virtual ~Solver();

    /** \brief Solve a standard (non-blocked) boundary integral equation.
      *
      * This function solves a boundary integral equation with given right-hand
      * side <tt>rhs</tt> of type <tt>GridFunction</tt> and returns a new <tt>Solution<tt>
      * object.
      *
      * \param[in] rhs
      * <tt>GridFunction</tt> representing the right-hand side function of the boundary
      * integral equation.
      *
      * \return A new <tt>Solution<tt> object, containing the solution of the boundary
      * integral equation.
      *
      */

    Solution<BasisFunctionType, ResultType> solve(
            const GridFunction<BasisFunctionType, ResultType>& rhs) const {
        return solveImplNonblocked(rhs); 
    }

    /** \brief Solve a block-operator system of boundary integral equations.
      *
      * This function solves a block system of boundary integral equations. It takes a
      * <tt>vector</tt> of variables of type <tt>GridFunction</tt> as its input.
      *
      * \param[in] rhs
      * <tt>vector</tt> of variables of type <tt>GridFunction</tt>
      *
      * \return A new <tt>BlockedSolution</tt> object, containing the solution of the system of
      * boundary integral equation.
      *
      */

    BlockedSolution<BasisFunctionType, ResultType> solve(
            const std::vector<GridFunction<BasisFunctionType, ResultType> >&
            rhs) const {
        return solveImplBlocked(rhs); 
    }

protected:
    static void checkConsistency(
        const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
        const GridFunction<BasisFunctionType, ResultType>& rhs,
        ConvergenceTestMode::Mode mode);
    static void checkConsistency(
        const BlockedBoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
        const std::vector<GridFunction<BasisFunctionType, ResultType> >& rhs,
        ConvergenceTestMode::Mode mode);

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
