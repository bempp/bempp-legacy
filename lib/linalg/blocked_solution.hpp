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

#ifndef bempp_blocked_solution_hpp
#define bempp_blocked_solution_hpp

#include "config_trilinos.hpp"

#include "solution_base.hpp"

#include "../assembly/grid_function.hpp"
#include <vector>

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType>
class BlockedSolution : public SolutionBase<BasisFunctionType, ResultType>
{
    typedef SolutionBase<BasisFunctionType, ResultType> Base;
public:
    typedef typename Base::MagnitudeType MagnitudeType;

    BlockedSolution(
            const std::vector<GridFunction<BasisFunctionType, ResultType> >& gridFunctions,
            const Thyra::SolveStatus<MagnitudeType> status);

    size_t gridFunctionCount() const;
    GridFunction<BasisFunctionType, ResultType>& gridFunction(size_t i);
    const GridFunction<BasisFunctionType, ResultType>& gridFunction(size_t i) const;

private:
    std::vector<GridFunction<BasisFunctionType, ResultType> > m_gridFunctions;
};

} // namespace Bempp

#endif
