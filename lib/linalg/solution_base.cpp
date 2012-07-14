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

#include "solution_base.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
SolutionBase<BasisFunctionType, ResultType>::SolutionBase(
    const Thyra::SolveStatus<MagnitudeType> status) :
    m_status(status)
{
}

template <typename BasisFunctionType, typename ResultType>
typename SolutionBase<BasisFunctionType, ResultType>::Status
SolutionBase<BasisFunctionType, ResultType>::status() const
{
    switch (m_status.solveStatus) {
    case Thyra::SOLVE_STATUS_CONVERGED:
        return CONVERGED;
    case Thyra::SOLVE_STATUS_UNCONVERGED:
        return UNCONVERGED;
    default:
        return UNKNOWN;
    }
}

template <typename BasisFunctionType, typename ResultType>
typename SolutionBase<BasisFunctionType, ResultType>::MagnitudeType
SolutionBase<BasisFunctionType, ResultType>::achievedTolerance() const
{
    return m_status.achievedTol;
}

template <typename BasisFunctionType, typename ResultType>
std::string 
SolutionBase<BasisFunctionType, ResultType>::solverMessage() const
{
    return m_status.message;
}

template <typename BasisFunctionType, typename ResultType>
Thyra::SolveStatus<typename SolutionBase<BasisFunctionType, ResultType>::MagnitudeType>
SolutionBase<BasisFunctionType, ResultType>::thyraSolveStatus() const
{
    return m_status;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(SolutionBase);

} // namespace Bempp
