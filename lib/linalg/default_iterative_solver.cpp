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

#include "default_iterative_solver.hpp"

#include "config_trilinos.hpp"

#ifdef WITH_TRILINOS

#include "Thyra_DefaultSpmdMultiVector.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_BelosLinearOpWithSolveFactory_decl.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerboseObject.hpp"


namespace Bempp {

template<typename ValueType>
DefaultIterativeSolver<ValueType>::DefaultIterativeSolver(DiscreteScalarValuedLinearOperator<ValueType>& discreteOperator,
                   DiscreteScalarValuedSourceTerm<ValueType>& rhs):
    m_discreteOperator(discreteOperator),
    m_rhs(&rhs,false){

    if (m_discreteOperator.range()->dim()!=rhs.range()->dim()) throw std::runtime_error("Dimension of LHS != Dimension of RHS.");

    }

template<typename ValueType>
DefaultIterativeSolver<ValueType>::DefaultIterativeSolver(DiscreteScalarValuedLinearOperator<ValueType>& discreteOperator,
                   std::vector<DiscreteScalarValuedSourceTerm<ValueType>* >& rhs) :
    m_discreteOperator(discreteOperator) {

    const size_t nrhs = rhs.size();

    if (nrhs<1) throw std::runtime_error("At least one right-hand side needed.");
    const size_t size = rhs[0]->range()->dim();

    // Check ranges

    for (size_t i=0;i<nrhs;i++)
        if (rhs[i]->range()->dim()!=size) throw std::runtime_error("Not all right-hand sides have same lengths.");

    if (m_discreteOperator.range()->dim()!=size) throw std::runtime_error("Dimension of LHS != Dimension of RHS.");

    // Create Multivector for RHS

    Teuchos::ArrayRCP<ValueType> data(nrhs*size);

    for (size_t i=0;i<nrhs;i++)
        for (size_t j=0;j<size;j++)
            data[i*nrhs+j]=rhs[i]->getPtr()[j];

    m_rhs=Teuchos::RCP<Thyra::MultiVectorBase<ValueType> >(new Thyra::DefaultSpmdMultiVector<ValueType>(
                                                   Thyra::defaultSpmdVectorSpace<ValueType>(size),
                                                   Thyra::defaultSpmdVectorSpace<ValueType>(nrhs),
                                                   data));

}

template<typename ValueType>
void DefaultIterativeSolver<ValueType>::addPreconditioner(Teuchos::RCP<Thyra::PreconditionerBase<ValueType> > preconditioner){
    m_preconditioner=preconditioner;
}

template<typename ValueType>
void DefaultIterativeSolver<ValueType>::initializeSolver(Teuchos::RCP<Teuchos::ParameterList> paramList){

    Teuchos::RCP<Teuchos::FancyOStream> out =Teuchos::VerboseObjectBase::getDefaultOStream();

    Thyra::BelosLinearOpWithSolveFactory<ValueType> invertibleOpFactory;
    invertibleOpFactory.setParameterList(paramList);
    invertibleOpFactory.setOStream(out);
    invertibleOpFactory.setVerbLevel(Teuchos::VERB_DEFAULT);
    // Wrap discreteOperator in reference counted Trilinos pointer.
    Teuchos::RCP<const Thyra::LinearOpBase<ValueType> > trilinosDiscreteLhs(&m_discreteOperator,false /*don't own*/);
    Teuchos::RCP<const Thyra::LinearOpSourceBase<ValueType> > linearOpSourcePtr(new Thyra::DefaultLinearOpSource<ValueType>(trilinosDiscreteLhs));

    if (!m_preconditioner.is_null()){
        // Preconditioner defined
        m_lhs = invertibleOpFactory.createOp();
        invertibleOpFactory.initializePreconditionedOp(linearOpSourcePtr,
                                                       m_preconditioner,
                                                       m_lhs.get(),
                                                       Thyra::SUPPORT_SOLVE_UNSPECIFIED);
        }
        else{
        // No preconditioner
        m_lhs=Thyra::linearOpWithSolve(invertibleOpFactory, trilinosDiscreteLhs);
    }


}

Teuchos::RCP<Teuchos::ParameterList> defaultGmresParameterList(double tol){

     Teuchos::RCP<Teuchos::ParameterList> paramList(new Teuchos::ParameterList("DefaultParameters"));
     paramList->set("Solver Type","Pseudo Block GMRES");
     Teuchos::ParameterList& solverTypesList = paramList->sublist("Solver Types");
     Teuchos::ParameterList& pseudoBlockGmresList = solverTypesList.sublist("Pseudo Block GMRES");
     pseudoBlockGmresList.set("Convergence Tolerance",tol);
     return paramList;

}

Teuchos::RCP<Teuchos::ParameterList> defaultCgParameterList(double tol){

     Teuchos::RCP<Teuchos::ParameterList> paramList(new Teuchos::ParameterList("DefaultParameters"));
     paramList->set("Solver Type","Pseudo Block CG");
     Teuchos::ParameterList& solverTypesList = paramList->sublist("Solver Types");
     Teuchos::ParameterList& pseudoBlockCgList = solverTypesList.sublist("Pseudo Block CG");
     pseudoBlockCgList.set("Convergence Tolerance",tol);
     return paramList;

}

#ifdef COMPILE_FOR_FLOAT
template class DefaultIterativeSolver<float>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class DefaultIterativeSolver<double>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class DefaultIterativeSolver<std::complex<float> >;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class DefaultIterativeSolver<std::complex<double> >;
#endif


}

#endif
