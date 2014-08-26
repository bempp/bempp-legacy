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

#include "bempp/common/config_trilinos.hpp"

#ifdef WITH_TRILINOS

#include "belos_solver_wrapper.hpp"

#include "real_wrapper_of_complex_thyra_linear_operator.hpp"
#include "real_wrapper_of_complex_thyra_preconditioner.hpp"
#include "../fiber/explicit_instantiation.hpp"

#include <Teuchos_ArrayRCP.hpp>
#include <Thyra_BelosLinearOpWithSolveFactory_decl.hpp>
#include <Thyra_DefaultLinearOpSource.hpp>
#include <Thyra_DetachedMultiVectorView.hpp>
#include <Thyra_LinearOpSourceBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_SolveSupportTypes.hpp>
#include <Thyra_OperatorVectorTypes.hpp>
#include <Thyra_PreconditionerBase.hpp>

#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

namespace Bempp {

// Overloaded template helper functions
namespace {

// Real ValueType
template <typename ValueType>
typename boost::enable_if<
    boost::is_same<ValueType, typename ScalarTraits<ValueType>::RealType>,
    Teuchos::RCP<Thyra::LinearOpWithSolveBase<
        typename ScalarTraits<ValueType>::RealType>>>::type
makeOperatorWithSolve(
    const Teuchos::RCP<Teuchos::ParameterList> &paramList,
    const Teuchos::RCP<const Thyra::LinearOpBase<ValueType>> &linOp,
    const Teuchos::RCP<const Thyra::PreconditionerBase<ValueType>> &
        preconditioner) {
  Teuchos::RCP<Teuchos::FancyOStream> out =
      Teuchos::VerboseObjectBase::getDefaultOStream();

  Thyra::BelosLinearOpWithSolveFactory<ValueType> invertibleOpFactory;
  invertibleOpFactory.setParameterList(paramList);
  invertibleOpFactory.setOStream(out);
  invertibleOpFactory.setVerbLevel(Teuchos::VERB_DEFAULT);

  Teuchos::RCP<Thyra::LinearOpWithSolveBase<ValueType>> result;
  if (preconditioner.is_null())
    // No preconditioner
    result = Thyra::linearOpWithSolve(invertibleOpFactory, linOp);
  else {
    // Preconditioner defined
    result = invertibleOpFactory.createOp();
    Teuchos::RCP<const Thyra::LinearOpSourceBase<ValueType>> linOpSourcePtr(
        new Thyra::DefaultLinearOpSource<ValueType>(linOp));
    invertibleOpFactory.initializePreconditionedOp(
        linOpSourcePtr, preconditioner, result.get(),
        Thyra::SUPPORT_SOLVE_UNSPECIFIED);
  }
  return result;
}

// Complex ValueType
template <typename ValueType>
typename boost::disable_if<
    boost::is_same<ValueType, typename ScalarTraits<ValueType>::RealType>,
    Teuchos::RCP<Thyra::LinearOpWithSolveBase<
        typename ScalarTraits<ValueType>::RealType>>>::type
makeOperatorWithSolve(
    const Teuchos::RCP<Teuchos::ParameterList> &paramList,
    const Teuchos::RCP<const Thyra::LinearOpBase<ValueType>> &linOp,
    const Teuchos::RCP<const Thyra::PreconditionerBase<ValueType>> &
        preconditioner) {
  typedef typename ScalarTraits<ValueType>::RealType RealType;

  Teuchos::RCP<Teuchos::FancyOStream> out =
      Teuchos::VerboseObjectBase::getDefaultOStream();

  Thyra::BelosLinearOpWithSolveFactory<RealType> invertibleOpFactory;
  invertibleOpFactory.setParameterList(paramList);
  invertibleOpFactory.setOStream(out);
  invertibleOpFactory.setVerbLevel(Teuchos::VERB_DEFAULT);

  Teuchos::RCP<const Thyra::LinearOpBase<RealType>> realLinOp(
      new RealWrapperOfComplexThyraLinearOperator<RealType>(linOp));

  Teuchos::RCP<Thyra::LinearOpWithSolveBase<RealType>> result;
  if (preconditioner.is_null())
    // No preconditioner
    result = Thyra::linearOpWithSolve(invertibleOpFactory, realLinOp);
  else {
    // Preconditioner defined
    result = invertibleOpFactory.createOp();
    Teuchos::RCP<const Thyra::LinearOpSourceBase<RealType>> realLinOpSourcePtr(
        new Thyra::DefaultLinearOpSource<RealType>(realLinOp));
    Teuchos::RCP<const Thyra::PreconditionerBase<RealType>> realPreconditioner(
        new RealWrapperOfComplexThyraPreconditioner<RealType>(preconditioner));
    invertibleOpFactory.initializePreconditionedOp(
        realLinOpSourcePtr, realPreconditioner, result.get(),
        Thyra::SUPPORT_SOLVE_UNSPECIFIED);
  }
  return result;
}

// Real ValueType
template <typename ValueType>
typename boost::enable_if<
    boost::is_same<ValueType, typename ScalarTraits<ValueType>::RealType>,
    Thyra::SolveStatus<typename ScalarTraits<ValueType>::RealType>>::type
reallySolve(const Thyra::LinearOpWithSolveBase<
                typename ScalarTraits<ValueType>::RealType> &op,
            const Thyra::EOpTransp trans,
            const Thyra::MultiVectorBase<ValueType> &rhs,
            const Teuchos::Ptr<Thyra::MultiVectorBase<ValueType>> &sol) {
  return op.solve(trans, rhs, sol);
}

// Complex ValueType
template <typename ValueType>
typename boost::disable_if<
    boost::is_same<ValueType, typename ScalarTraits<ValueType>::RealType>,
    Thyra::SolveStatus<typename ScalarTraits<ValueType>::RealType>>::type
reallySolve(const Thyra::LinearOpWithSolveBase<
                typename ScalarTraits<ValueType>::RealType> &op,
            const Thyra::EOpTransp trans,
            const Thyra::MultiVectorBase<ValueType> &rhs,
            const Teuchos::Ptr<Thyra::MultiVectorBase<ValueType>> &sol) {
  typedef typename ScalarTraits<ValueType>::RealType MagnitudeType;

  // Wrap the right-hand-side data in a real multivector of twice as many rows

  const Teuchos::Ordinal rhsRowCount = rhs.range()->dim();
  const Teuchos::Ordinal rhsColCount = rhs.domain()->dim();

  Teuchos::ArrayRCP<const MagnitudeType> realRhsArray;
  RTOpPack::ConstSubMultiVectorView<MagnitudeType> realRhsView;
  Teuchos::RCP<const Thyra::MultiVectorBase<MagnitudeType>> realRhs;

  Thyra::ConstDetachedMultiVectorView<ValueType> rhsView(
      Teuchos::rcpFromRef(rhs));
  if (rhsView.leadingDim() == rhsRowCount) {
    // contiguous vector
    const MagnitudeType *realRhsData =
        reinterpret_cast<const MagnitudeType *>(rhsView.values());
    realRhsArray =
        Teuchos::arcp(realRhsData, 0,                        // lowerOffset
                      2 * rhsRowCount * rhsColCount, false); // doesn't own
    realRhsView.initialize(0,                                // globalOffset
                           2 * rhsRowCount, 0,               // colOffset
                           rhsColCount, realRhsArray,
                           2 * rhsRowCount); // leadingDim
    realRhs = Thyra::createMembersView(op.domain(), realRhsView);
  } else
    throw std::runtime_error("reallySolve(): discontiguous multivectors "
                             "are not supported yet");

  // Wrap the solution data in a real multivector of twice as many rows

  const Teuchos::Ordinal solRowCount = sol->range()->dim();
  const Teuchos::Ordinal solColCount = sol->domain()->dim();

  Teuchos::ArrayRCP<MagnitudeType> realSolArray;
  RTOpPack::SubMultiVectorView<MagnitudeType> realSolView;
  Teuchos::RCP<Thyra::MultiVectorBase<MagnitudeType>> realSol;

  Thyra::DetachedMultiVectorView<ValueType> solView(Teuchos::rcpFromPtr(sol));
  if (solView.leadingDim() == solRowCount) {
    // contiguous vector
    MagnitudeType *realSolData =
        reinterpret_cast<MagnitudeType *>(solView.values());
    realSolArray =
        Teuchos::arcp(realSolData, 0,                        // lowerOffset
                      2 * solRowCount * solColCount, false); // doesn't own
    realSolView.initialize(0,                                // globalOffset
                           2 * solRowCount, 0,               // colOffset
                           solColCount, realSolArray,
                           2 * solRowCount); // leadingDim
    realSol = Thyra::createMembersView(op.range(), realSolView);
  } else
    throw std::runtime_error("reallySolve(): discontiguous multivectors "
                             "are not supported yet");

  return op.solve(trans, *realRhs, realSol.ptr());
}

} // namespace

// BelosSolverWrapper member functions

template <typename ValueType>
BelosSolverWrapper<ValueType>::BelosSolverWrapper(
    const Teuchos::RCP<const Thyra::LinearOpBase<ValueType>> &linOp)
    : m_linOp(linOp) {}

template <typename ValueType>
BelosSolverWrapper<ValueType>::~BelosSolverWrapper() {}

template <typename ValueType>
void BelosSolverWrapper<ValueType>::setPreconditioner(const Teuchos::RCP<
    const Thyra::PreconditionerBase<ValueType>> &preconditioner) {
  m_preconditioner = preconditioner;
}

template <typename ValueType>
void BelosSolverWrapper<ValueType>::initializeSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &paramList) {
  m_linOpWithSolve =
      makeOperatorWithSolve(paramList, m_linOp, m_preconditioner);
}

template <typename ValueType>
Thyra::SolveStatus<typename BelosSolverWrapper<ValueType>::MagnitudeType>
BelosSolverWrapper<ValueType>::solve(
    const Thyra::EOpTransp trans, const Thyra::MultiVectorBase<ValueType> &rhs,
    const Teuchos::Ptr<Thyra::MultiVectorBase<ValueType>> &sol) const {
  if (m_linOpWithSolve.is_null())
    throw std::runtime_error("BelosSolverWrapper::solve(): "
                             "solver not initialized");
  return reallySolve(*m_linOpWithSolve, trans, rhs, sol);
}

// nonmember functions

namespace {

template <typename MagnitudeType>
Teuchos::RCP<Teuchos::ParameterList> inline defaultGmresParameterListInternal(
    MagnitudeType tol, int maxIterationCount) {
  Teuchos::RCP<Teuchos::ParameterList> paramList(
      new Teuchos::ParameterList("DefaultParameters"));
  paramList->set("Solver Type", "Pseudo Block GMRES");
  Teuchos::ParameterList &solverTypesList = paramList->sublist("Solver Types");
  Teuchos::ParameterList &pseudoBlockGmresList =
      solverTypesList.sublist("Pseudo Block GMRES");
  pseudoBlockGmresList.set("Convergence Tolerance", tol);
  pseudoBlockGmresList.set("Maximum Iterations", maxIterationCount);
  return paramList;
}

template <typename MagnitudeType>
Teuchos::RCP<Teuchos::ParameterList> inline defaultCgParameterListInternal(
    MagnitudeType tol, int maxIterationCount) {
  Teuchos::RCP<Teuchos::ParameterList> paramList(
      new Teuchos::ParameterList("DefaultParameters"));
  paramList->set("Solver Type", "Pseudo Block CG");
  Teuchos::ParameterList &solverTypesList = paramList->sublist("Solver Types");
  Teuchos::ParameterList &pseudoBlockCgList =
      solverTypesList.sublist("Pseudo Block CG");
  pseudoBlockCgList.set("Convergence Tolerance", tol);
  pseudoBlockCgList.set("Maximum Iterations", maxIterationCount);
  return paramList;
}

} // namespace

Teuchos::RCP<Teuchos::ParameterList>
defaultGmresParameterList(double tol, int maxIterationCount) {
  return defaultGmresParameterListInternal(tol, maxIterationCount);
}

Teuchos::RCP<Teuchos::ParameterList>
defaultGmresParameterList(float tol, int maxIterationCount) {
  return defaultGmresParameterListInternal(tol, maxIterationCount);
}

Teuchos::RCP<Teuchos::ParameterList>
defaultCgParameterList(double tol, int maxIterationCount) {
  return defaultCgParameterListInternal(tol, maxIterationCount);
}

Teuchos::RCP<Teuchos::ParameterList>
defaultCgParameterList(float tol, int maxIterationCount) {
  return defaultCgParameterListInternal(tol, maxIterationCount);
}

FIBER_INSTANTIATE_CLASS_TEMPLATED_ON_RESULT(BelosSolverWrapper);

} // namespace Bempp

#endif // WITH_TRILINOS
