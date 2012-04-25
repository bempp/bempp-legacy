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

#include "default_direct_solver.hpp"

#include "../assembly/linear_operator.hpp"
#include "../assembly/discrete_linear_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

template <typename ArgumentType, typename ResultType>
DefaultDirectSolver<ArgumentType, ResultType>::DefaultDirectSolver(
        const LinearOperator<ArgumentType, ResultType>& linearOperator,
        const GridFunction<ArgumentType, ResultType>& gridFunction) :
    m_linearOperator(linearOperator), m_gridFunction(gridFunction),
    m_solution(), m_status(Solver<ArgumentType, ResultType>::UNKNOWN)
{
    if (!linearOperator.isAssembled())
        throw std::runtime_error("DefaultDirectSolver::DefaultDirectSolver(): "
                                 "operator is not assembled");
    if (&linearOperator.trialSpace() != &gridFunction.space())
        throw std::runtime_error("DefaultDirectSolver::DefaultDirectSolver(): "
                                 "spaces do not match");
}

template <typename ArgumentType, typename ResultType>
void DefaultDirectSolver<ArgumentType, ResultType>::solve()
{
    m_solution = arma::solve(
                m_linearOperator.assembledDiscreteLinearOperator().asMatrix(),
                m_gridFunction.coefficients().asArmadilloVector());
    m_status = Solver<ArgumentType, ResultType>::CONVERGED;
}

template <typename ArgumentType, typename ResultType>
GridFunction<ArgumentType, ResultType> DefaultDirectSolver<ArgumentType, ResultType>::getResult() const
{
    return GridFunction<ArgumentType, ResultType>(
                m_linearOperator.testSpace(),
                m_solution);
}

template <typename ArgumentType, typename ResultType>
typename Solver<ArgumentType, ResultType>::EStatus DefaultDirectSolver<ArgumentType, ResultType>::getStatus() const
{
    return m_status;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_KERNEL(DefaultDirectSolver);

} // namespace Bempp
