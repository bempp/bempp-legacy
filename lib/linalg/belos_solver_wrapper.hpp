#ifndef belos_solver_wrapper_hpp
#define belos_solver_wrapper_hpp

#include "../common/common.hpp"

#include "config_trilinos.hpp"

#ifdef WITH_TRILINOS

#include "../common/scalar_traits.hpp"

#include <Thyra_SolveSupportTypes.hpp>
#include <Thyra_PreconditionerBase.hpp>

namespace Bempp
{

template <typename ValueType>
class BelosSolverWrapper
{
public:
    typedef typename ScalarTraits<ValueType>::RealType MagnitudeType;

    BelosSolverWrapper(
            const Teuchos::RCP<const Thyra::LinearOpBase<ValueType> >& linOp);

    void addPreconditioner(
            const Teuchos::RCP<const Thyra::PreconditionerBase<ValueType> >& preconditioner);

    void initializeSolver(
            const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    Thyra::SolveStatus<MagnitudeType> solve(
            const Thyra::EOpTransp trans,
            const Thyra::MultiVectorBase<ValueType>& rhs,
            const Teuchos::Ptr<Thyra::MultiVectorBase<ValueType> >& sol) const;

private:
    Teuchos::RCP<const Thyra::LinearOpBase<ValueType> > m_linOp;
    Teuchos::RCP<const Thyra::PreconditionerBase<ValueType> > m_preconditioner;
    Teuchos::RCP<const Thyra::LinearOpWithSolveBase<MagnitudeType> > m_linOpWithSolve;
};

Teuchos::RCP<Teuchos::ParameterList> defaultGmresParameterList(double tol);
Teuchos::RCP<Teuchos::ParameterList> defaultCgParameterList(double tol);

} // namespace Bempp

#endif // WITH_TRILINOS

#endif
