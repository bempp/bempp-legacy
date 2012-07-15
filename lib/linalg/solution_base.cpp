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

#include "config_trilinos.hpp"

#include "solution_base.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp
{

#ifdef WITH_TRILINOS
template <typename BasisFunctionType, typename ResultType>
SolutionBase<BasisFunctionType, ResultType>::SolutionBase(
        const Thyra::SolveStatus<MagnitudeType> status) :
    m_achievedTolerance(status.achievedTol),
    m_message(status.message),
    m_extraParameters(status.extraParameters)
{
    switch (status.solveStatus) {
    case Thyra::SOLVE_STATUS_CONVERGED:
        m_status = CONVERGED; break;
    case Thyra::SOLVE_STATUS_UNCONVERGED:
        m_status = UNCONVERGED; break;
    default:
        m_status = UNKNOWN;
    }
}
#endif // WITH_TRILINOS

template <typename BasisFunctionType, typename ResultType>
SolutionBase<BasisFunctionType, ResultType>::SolutionBase(
        Status status, MagnitudeType achievedTolerance, std::string message) :
    m_status(status), 
    m_achievedTolerance(achievedTolerance),
    m_message(message)
{
}

template <typename BasisFunctionType, typename ResultType>
typename SolutionBase<BasisFunctionType, ResultType>::Status
SolutionBase<BasisFunctionType, ResultType>::status() const
{
    return m_status;
}

template <typename BasisFunctionType, typename ResultType>
typename SolutionBase<BasisFunctionType, ResultType>::MagnitudeType
SolutionBase<BasisFunctionType, ResultType>::achievedTolerance() const
{
    return m_achievedTolerance;
}

template <typename BasisFunctionType, typename ResultType>
std::string 
SolutionBase<BasisFunctionType, ResultType>::solverMessage() const
{
    return m_message;
}

#ifdef WITH_TRILINOS
template <typename BasisFunctionType, typename ResultType>
Thyra::RCP<Teuchos::ParameterList> 
SolutionBase<BasisFunctionType, ResultType>::extraParameters() const
{
    return m_extraParameters;
}
#endif // WITH_TRILINOS

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(SolutionBase);

} // namespace Bempp
