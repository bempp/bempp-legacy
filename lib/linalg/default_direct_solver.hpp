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

namespace Bempp
{

template <typename ValueType> class LinearOperator;

template <typename ValueType>
class DefaultDirectSolver : public Solver<ValueType>
{
public:
    DefaultDirectSolver(const LinearOperator<ValueType>& linOp,
                        const GridFunction<ValueType>& gridFun);

    virtual void solve();

    virtual GridFunction<ValueType> getResult() const;
    virtual typename Solver<ValueType>::EStatus getStatus() const;

private:
    const LinearOperator<ValueType>& m_linearOperator;
    const GridFunction<ValueType>& m_gridFunction;
    arma::Col<ValueType> m_solution;
    typename Solver<ValueType>::EStatus m_status;
};

} // namespace Bempp

#endif
