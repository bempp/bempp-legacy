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

#ifndef bempp_default_iterative_solver_hpp
#define bempp_default_iterative_solver_hpp

#include "../common/common.hpp"

#include "config_trilinos.hpp"

#ifdef WITH_TRILINOS

// #include "solver.hpp"

#include "solution.hpp"
#include "belos_solver_wrapper_fwd.hpp" // for default parameter lists
#include "blocked_solution.hpp"
#include "../common/armadillo_fwd.hpp"
#include "../common/scalar_traits.hpp"
#include "../common/shared_ptr.hpp"

#include <boost/scoped_ptr.hpp>

namespace Thyra
{
template <typename ValueType> class PreconditionerBase;
template <typename ValueType> class SolveStatus;
}

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType> class BoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class BlockedBoundaryOperator;
template <typename BasisFunctionType, typename ResultType> class GridFunction;

template <typename BasisFunctionType, typename ResultType>
class DefaultIterativeSolver // : public Solver<BasisFunctionType, ResultType>
{
public:
    typedef typename ScalarTraits<ResultType>::RealType MagnitudeType;

    enum ConvergenceTestMode {
        TEST_CONVERGENCE_IN_DUAL_TO_RANGE,
        TEST_CONVERGENCE_IN_RANGE
    };

    DefaultIterativeSolver(
        const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
            ConvergenceTestMode mode = TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
    DefaultIterativeSolver(
        const BlockedBoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
        ConvergenceTestMode mode = TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
    virtual ~DefaultIterativeSolver();

    void setPreconditioner(
            const Teuchos::RCP<const Thyra::PreconditionerBase<ResultType> >& preconditioner);
    void initializeSolver(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    virtual Solution<BasisFunctionType, ResultType> solve(
            const GridFunction<BasisFunctionType, ResultType>& rhs) const;
    virtual BlockedSolution<BasisFunctionType, ResultType> solve(
            const std::vector<GridFunction<BasisFunctionType, ResultType> >&
            rhs) const;

private:
    struct Impl;
    boost::scoped_ptr<Impl> m_impl;
};

} // namespace Bempp

#endif // WITH_TRILINOS

#endif
