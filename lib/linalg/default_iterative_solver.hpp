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

#ifndef bempp_default_iterative_solver_hpp
#define bempp_default_iterative_solver_hpp

#include "../common/common.hpp"

#include "config_trilinos.hpp"

#ifdef WITH_TRILINOS

#include "solver.hpp"

#include "belos_solver_wrapper.hpp"

#include "../common/armadillo_fwd.hpp"
#include "../common/shared_ptr.hpp"

namespace Bempp
{

template <typename BasisFunctionType, typename ResultType> class BoundaryOperator;

template <typename BasisFunctionType, typename ResultType>
class DefaultIterativeSolver : public Solver<BasisFunctionType, ResultType>
{
public:
    typedef typename ScalarTraits<ResultType>::RealType MagnitudeType;

    DefaultIterativeSolver(const BoundaryOperator<BasisFunctionType, ResultType>& boundaryOp,
                           const GridFunction<BasisFunctionType, ResultType>& rhsGridFun);

    void addPreconditioner(
            const Teuchos::RCP<const Thyra::PreconditionerBase<ResultType> >& preconditioner);
    void initializeSolver(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    virtual void solve();

    virtual GridFunction<BasisFunctionType, ResultType> getResult() const;
    virtual typename Solver<BasisFunctionType, ResultType>::EStatus getStatus() const;
    MagnitudeType getSolveTolerance() const;
    std::string getSolverMessage() const;
    Thyra::SolveStatus<MagnitudeType> getThyraSolveStatus() const;

private:
    shared_ptr<const Context<BasisFunctionType, ResultType> > m_context;
    BelosSolverWrapper<ResultType> m_belosSolverWrapper;
    shared_ptr<const Space<BasisFunctionType> > m_space;
    Teuchos::RCP<Thyra::MultiVectorBase<ResultType> > m_rhs;

    // as long as a GridFunction is initialised with an Armadillo array
    // rather than a Vector, this is the right format for solution storage
    arma::Col<ResultType> m_sol;
    Thyra::SolveStatus<MagnitudeType> m_status;
};

} // namespace Bempp

#endif // WITH_TRILINOS

#endif
