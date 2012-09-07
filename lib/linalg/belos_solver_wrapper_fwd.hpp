#ifndef bempp_belos_solver_wrapper_fwd_hpp
#define bempp_belos_solver_wrapper_fwd_hpp

#include "../common/common.hpp"

#include "bempp/common/config_trilinos.hpp"

#ifdef WITH_TRILINOS
#include <Thyra_SolveSupportTypes.hpp>

namespace Bempp
{

/** \cond FORWARD_DECL */
template <typename ValueType> class BelosSolverWrapper;
/** \endcond */

Teuchos::RCP<Teuchos::ParameterList> defaultGmresParameterList(
        double tol, int maxIterationCount = 1000);
Teuchos::RCP<Teuchos::ParameterList> defaultCgParameterList(
        double tol, int maxIterationCount = 1000);
Teuchos::RCP<Teuchos::ParameterList> defaultGmresParameterList(
        float tol, int maxIterationCount = 1000);
Teuchos::RCP<Teuchos::ParameterList> defaultCgParameterList(
        float tol, int maxIterationCount = 1000);

} // namespace Bempp

#endif // WITH_TRILINOS

#endif
