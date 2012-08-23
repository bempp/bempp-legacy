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

#include "bempp/common/config_trilinos.hpp"

#include "solution_base.hpp"

#include "../assembly/grid_function.hpp"
#include <vector>

namespace Bempp
{

/** \ingroup linalg
  * \brief This class holds the solution of a block operator system
           together with various information about the solution.
  *
  */

template <typename BasisFunctionType, typename ResultType>
class BlockedSolution : public SolutionBase<BasisFunctionType, ResultType>
{
public:
    typedef SolutionBase<BasisFunctionType, ResultType> Base;
    typedef typename Base::MagnitudeType MagnitudeType;

#ifdef WITH_TRILINOS
    /** \brief Constructor */
    BlockedSolution(
            const std::vector<GridFunction<BasisFunctionType, ResultType> >& gridFunctions,
            const Thyra::SolveStatus<MagnitudeType> status);
#endif // WITH_TRILINOS
    /** \brief Constructor */
    BlockedSolution(
            const std::vector<GridFunction<BasisFunctionType, ResultType> >& gridFunctions,
            SolutionStatus::Status status,
            MagnitudeType achievedTolerance = Base::unknownTolerance(),
            std::string message = "");

    size_t gridFunctionCount() const;
    GridFunction<BasisFunctionType, ResultType>& gridFunction(size_t i);
    const GridFunction<BasisFunctionType, ResultType>& gridFunction(size_t i) const;

private:
    std::vector<GridFunction<BasisFunctionType, ResultType> > m_gridFunctions;
};

} // namespace Bempp

#endif
