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

#include "solver.hpp"

#include "config_trilinos.hpp"

#ifdef WITH_TRILINOS

#include <vector>
#include <armadillo>
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

#include "../assembly/discrete_linear_operator.hpp"
#include "../assembly/discrete_aca_linear_operator.hpp"
#include "../assembly/vector.hpp"

namespace Bempp
{

template <typename ValueType> class LinearOperator;

template <typename ValueType>
class DefaultIterativeSolver : public Solver<ValueType>
{
public:
    DefaultIterativeSolver(const LinearOperator<ValueType>& linOp,
                           const GridFunction<ValueType>& gridFun);

    void addPreconditioner(
            Teuchos::RCP<const Thyra::PreconditionerBase<ValueType> > preconditioner);
    void initializeSolver(Teuchos::RCP<Teuchos::ParameterList> paramList);

    virtual void solve();

    virtual GridFunction<ValueType> getResult() const;
    virtual typename Solver<ValueType>::EStatus getStatus() const;
    double getSolveTolerance() const;
    std::string getSolverMessage() const;
    Thyra::SolveStatus<ValueType> getThyraSolveStatus() const;

private:
    const DiscreteLinearOperator<ValueType>& m_discreteOperator;
    const Space<ValueType>& m_space;

    Teuchos::RCP<Thyra::LinearOpWithSolveBase<ValueType> > m_lhs;
    Teuchos::RCP<Thyra::MultiVectorBase<ValueType> > m_rhs;
    Teuchos::RCP<Thyra::MultiVectorBase<ValueType> > m_sol;
    Thyra::SolveStatus<ValueType> m_status;
    Teuchos::RCP<const Thyra::PreconditionerBase<ValueType> > m_preconditioner;
};

Teuchos::RCP<Teuchos::ParameterList> defaultGmresParameterList(double tol);
Teuchos::RCP<Teuchos::ParameterList> defaultCgParameterList(double tol);

} // namespace Bempp

#endif // WITH_TRILINOS
#endif
