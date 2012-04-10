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

#ifndef default_gmres_solver_hpp
#define default_gmres_solver_hpp

#include "config_trilinos.hpp"

#ifdef WITH_TRILINOS

#include <vector>
#include <armadillo>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Thyra_PreconditionerBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include <Thyra_BelosLinearOpWithSolveFactory_decl.hpp>
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_PreconditionerBase.hpp"

#include "../assembly/discrete_scalar_valued_linear_operator.hpp"
#include "../assembly/discrete_aca_scalar_valued_linear_operator.hpp"
#include "../assembly/discrete_scalar_valued_source_term.hpp"



namespace Bempp {

enum EStatus {CONVERGED, UNCONVGERGED, UNKNOWN};

template <typename ValueType>
class DefaultIterativeSolver
{
public:


    DefaultIterativeSolver(DiscreteScalarValuedLinearOperator<ValueType>& discreteOperator,
                       DiscreteScalarValuedSourceTerm<ValueType>& rhs);
    DefaultIterativeSolver(DiscreteScalarValuedLinearOperator<ValueType>& discreteOperator,
                       std::vector<DiscreteScalarValuedSourceTerm<ValueType>* >& rhs);


    void addPreconditioner(Teuchos::RCP<const Thyra::PreconditionerBase<ValueType> > preconditioner);


    void initializeSolver(Teuchos::RCP<Teuchos::ParameterList> paramList);



    void solve();



    arma::Mat<ValueType> getResult();

    EStatus getStatus();

    double getSolveTol();

    std::string getSolverMessage();

    Thyra::SolveStatus<ValueType> getThyraSolveStatus();

private:

    DiscreteScalarValuedLinearOperator<ValueType>& m_discreteOperator;

    Teuchos::RCP<Thyra::LinearOpWithSolveBase<ValueType> > m_lhs;
    Teuchos::RCP<Thyra::MultiVectorBase<ValueType> > m_rhs;
    Teuchos::RCP<Thyra::MultiVectorBase<ValueType> > m_sol;
    Thyra::SolveStatus<ValueType> m_status;
    Teuchos::RCP<const Thyra::PreconditionerBase<ValueType> > m_preconditioner;


};

Teuchos::RCP<Teuchos::ParameterList> defaultGmresParameterList(double tol);
Teuchos::RCP<Teuchos::ParameterList> defaultCgParameterList(double tol);

}

#endif // WITH_TRILINOS
#endif // default_gmres_solver_hpp
