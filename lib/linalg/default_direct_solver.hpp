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

#ifndef bempp_default_direct_solver_hpp
#define bempp_default_direct_solver_hpp

#include "solver.hpp"

#include <boost/scoped_ptr.hpp>

namespace Bempp
{

/** \ingroup linalg
  * \brief Default Interface to the solution of boundary integral equations using a dense LU decomposition.
  *
  * This class provides an interface to the direct solution of boundary integral equations using standard LU.
  * It is implemented via a call to the <a href="http://arma.sourceforge.net">Armadillo</a> dense solver, which
  * is an interface to the corresponding Lapack routines.
  *
  */

template <typename BasisFunctionType, typename ResultType>
class DefaultDirectSolver : public Solver<BasisFunctionType, ResultType>
{
public:
    typedef Solver<BasisFunctionType, ResultType> Base;

    /** \brief Construct a solver for a non-blocked boundary operator. */
    DefaultDirectSolver(
            const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOp);
    /** \brief Construct a solver for a blocked boundary operator. */
    DefaultDirectSolver(
            const BlockedBoundaryOperator<BasisFunctionType, ResultType>& boundaryOp);
    ~DefaultDirectSolver();

private:
    virtual Solution<BasisFunctionType, ResultType> solveImplNonblocked(
            const GridFunction<BasisFunctionType, ResultType>& rhs) const;
    virtual BlockedSolution<BasisFunctionType, ResultType> solveImplBlocked(
            const std::vector<GridFunction<BasisFunctionType, ResultType> >&
            rhs) const;

private:
    struct Impl;
    boost::scoped_ptr<Impl> m_impl;
};

} // namespace Bempp

#endif
