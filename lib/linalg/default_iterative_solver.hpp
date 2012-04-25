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

#include "config_trilinos.hpp"

#ifdef WITH_TRILINOS

#include "solver.hpp"

#include "../assembly/discrete_linear_operator.hpp"
#include "../assembly/discrete_aca_linear_operator.hpp"
#include "../assembly/vector.hpp"

#include <armadillo>
#include <vector>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Thyra_PreconditionerBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_BelosLinearOpWithSolveFactory_decl.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_MultiVectorBase.hpp>
#include <Thyra_PreconditionerBase.hpp>

namespace Bempp
{

template <typename ArgumentType, typename ResultType> class LinearOperator;

template <typename ArgumentType, typename ResultType>
class DefaultIterativeSolver : public Solver<ArgumentType, ResultType>
{
public:
    DefaultIterativeSolver(const LinearOperator<ArgumentType, ResultType>& linOp,
                           const GridFunction<ArgumentType, ResultType>& gridFun);

    void addPreconditioner(
            Teuchos::RCP<const Thyra::PreconditionerBase<ResultType> > preconditioner);
    void initializeSolver(Teuchos::RCP<Teuchos::ParameterList> paramList);

    virtual void solve();

    virtual GridFunction<ArgumentType, ResultType> getResult() const;
    virtual typename Solver<ArgumentType, ResultType>::EStatus getStatus() const;
    double getSolveTolerance() const;
    std::string getSolverMessage() const;
    Thyra::SolveStatus<ResultType> getThyraSolveStatus() const;

private:
    const DiscreteLinearOperator<ResultType>& m_discreteOperator;
    const Space<ArgumentType>& m_space;

    Teuchos::RCP<Thyra::LinearOpWithSolveBase<ResultType> > m_lhs;
    Teuchos::RCP<Thyra::MultiVectorBase<ResultType> > m_rhs;
    Teuchos::RCP<Thyra::MultiVectorBase<ResultType> > m_sol;
    Thyra::SolveStatus<ResultType> m_status;
    Teuchos::RCP<const Thyra::PreconditionerBase<ResultType> > m_preconditioner;
};

Teuchos::RCP<Teuchos::ParameterList> defaultGmresParameterList(double tol);
Teuchos::RCP<Teuchos::ParameterList> defaultCgParameterList(double tol);

} // namespace Bempp

#endif // WITH_TRILINOS
#endif
