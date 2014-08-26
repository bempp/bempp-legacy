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

#include "solution.hpp"
#include "../fiber/explicit_instantiation.hpp"

namespace Bempp {

template <typename BasisFunctionType, typename ResultType>
Solution<BasisFunctionType, ResultType>::Solution()
    : Base(SolutionStatus::UNKNOWN) {}

#ifdef WITH_TRILINOS
template <typename BasisFunctionType, typename ResultType>
Solution<BasisFunctionType, ResultType>::Solution(
    const GridFunction<BasisFunctionType, ResultType> &gridFunction,
    const Thyra::SolveStatus<MagnitudeType> status)
    : Base(status), m_gridFunction(gridFunction) {}
#endif // WITH_TRILINOS

template <typename BasisFunctionType, typename ResultType>
Solution<BasisFunctionType, ResultType>::Solution(
    const GridFunction<BasisFunctionType, ResultType> &gridFunction,
    SolutionStatus::Status status, MagnitudeType achievedTolerance,
    std::string message)
    : Base(status, achievedTolerance, message), m_gridFunction(gridFunction) {}

template <typename BasisFunctionType, typename ResultType>
GridFunction<BasisFunctionType, ResultType> &
Solution<BasisFunctionType, ResultType>::gridFunction() {
  return m_gridFunction;
}

template <typename BasisFunctionType, typename ResultType>
const GridFunction<BasisFunctionType, ResultType> &
Solution<BasisFunctionType, ResultType>::gridFunction() const {
  return m_gridFunction;
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_BASIS_AND_RESULT(Solution);

} // namespace Bempp
