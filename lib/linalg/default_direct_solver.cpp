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

#include "../assembly/abstract_boundary_operator.hpp"
#include "../assembly/boundary_operator.hpp"
#include "../assembly/discrete_boundary_operator.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
DefaultDirectSolver<BasisFunctionType, ResultType>::DefaultDirectSolver(
        const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOperator,
        const GridFunction<BasisFunctionType, ResultType>& gridFunction) :
    m_boundaryOperator(boundaryOperator), m_gridFunction(gridFunction),
    m_solution(), m_status(Solver<BasisFunctionType, ResultType>::UNKNOWN)
{
    if (&boundaryOperator.abstractOperator()->domain() != &gridFunction.space())
        throw std::runtime_error("DefaultDirectSolver::DefaultDirectSolver(): "
                                 "spaces do not match");
}

template <typename BasisFunctionType, typename ResultType>
void DefaultDirectSolver<BasisFunctionType, ResultType>::solve()
{
    m_solution = arma::solve(
                m_boundaryOperator.weakForm()->asMatrix(),
                m_gridFunction.projections());
    m_status = Solver<BasisFunctionType, ResultType>::CONVERGED;
}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType>
DefaultDirectSolver<BasisFunctionType, ResultType>::getResult() const
{
    return gridFunctionFromCoefficients(
                m_boundaryOperator.context(),
                m_boundaryOperator.abstractOperator()->domain(),
                m_boundaryOperator.abstractOperator()->domain(), // is this the right choice?
                m_solution);
}

template <typename BasisFunctionType, typename ResultType>
typename Solver<BasisFunctionType, ResultType>::EStatus
DefaultDirectSolver<BasisFunctionType, ResultType>::getStatus() const
{
    return m_status;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(DefaultDirectSolver);

} // namespace Bempp
