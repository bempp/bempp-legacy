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

namespace Bempp
{

template<typename ValueType>
DefaultDirectSolver<ValueType>::DefaultDirectSolver(
        const LinearOperator<ValueType>& linearOperator,
        const GridFunction<ValueType>& gridFunction) :
    m_linearOperator(linearOperator), m_gridFunction(gridFunction),
    m_solution(), m_status(Solver<ValueType>::UNKNOWN)
{
    if (!linearOperator.isAssembled())
        throw std::runtime_error("DefaultDirectSolver::DefaultDirectSolver(): "
                                 "operator is not assembled");
    if (&linearOperator.trialSpace() != &gridFunction.space())
        throw std::runtime_error("DefaultDirectSolver::DefaultDirectSolver(): "
                                 "spaces do not match");
}

template <typename ValueType>
void DefaultDirectSolver<ValueType>::solve()
{
    m_solution = arma::solve(
                m_linearOperator.assembledDiscreteLinearOperator().asMatrix(),
                m_gridFunction.coefficients().asArmadilloVector());
    m_status = Solver<ValueType>::CONVERGED;
}

template <typename ValueType>
GridFunction<ValueType> DefaultDirectSolver<ValueType>::getResult() const
{
    return GridFunction<ValueType>(
                m_linearOperator.testSpace(),
                m_solution);
}

template <typename ValueType>
typename Solver<ValueType>::EStatus DefaultDirectSolver<ValueType>::getStatus() const
{
    return m_status;
}

#ifdef COMPILE_FOR_FLOAT
template class DefaultDirectSolver<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class DefaultDirectSolver<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class DefaultDirectSolver<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class DefaultDirectSolver<std::complex<double> >;
#endif

} // namespace Bempp
