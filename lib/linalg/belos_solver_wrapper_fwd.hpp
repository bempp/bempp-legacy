#ifndef bempp_belos_solver_wrapper_fwd_hpp
#define bempp_belos_solver_wrapper_fwd_hpp

#include "../common/common.hpp"

#include "config_trilinos.hpp"

#ifdef WITH_TRILINOS
#include <Thyra_SolveSupportTypes.hpp>

namespace Bempp
{

template <typename ValueType> class BelosSolverWrapper;

Teuchos::RCP<Teuchos::ParameterList> defaultGmresParameterList(double tol);
Teuchos::RCP<Teuchos::ParameterList> defaultCgParameterList(double tol);
Teuchos::RCP<Teuchos::ParameterList> defaultGmresParameterList(float tol);
Teuchos::RCP<Teuchos::ParameterList> defaultCgParameterList(float tol);

} // namespace Bempp

#endif // WITH_TRILINOS

#endif
