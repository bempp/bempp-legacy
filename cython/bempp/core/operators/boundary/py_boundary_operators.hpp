#ifndef bempp_py_boundary_operators_hpp
#define bempp_py_boundary_operators_hpp

#include "bempp/common/boost_make_shared_fwd.hpp"
#include "bempp/common/shared_ptr.hpp"
#include "bempp/assembly/general_elementary_local_operator_imp.hpp"
#include "bempp/assembly/abstract_boundary_operator.hpp"
#include "bempp/assembly/elementary_integral_operator.hpp"
#include "bempp/assembly/elementary_local_operator.hpp"
#include "bempp/assembly/symmetry.hpp"
#include "bempp/assembly/identity_operator.hpp"
#include "bempp/assembly/laplace_beltrami_3d_operator.hpp"
#include "bempp/assembly/maxwell_3d_identity_operator.hpp"
#include "bempp/assembly/maxwell_3d_single_layer_boundary_operator.hpp"
#include "bempp/assembly/maxwell_3d_double_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_double_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_single_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_adjoint_double_layer_boundary_operator.hpp"
#include "bempp/assembly/laplace_3d_hypersingular_boundary_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_double_layer_boundary_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_adjoint_double_layer_boundary_operator.hpp"
#include "bempp/assembly/modified_helmholtz_3d_hypersingular_boundary_operator.hpp"
#include "bempp/fiber/hdiv_function_value_functor.hpp"
#include "bempp/fiber/surface_curl_3d_functor.hpp"
#include "bempp/fiber/scalar_function_value_times_normal_functor.hpp"
#include "bempp/fiber/scalar_function_value_functor.hpp"
#include "bempp/fiber/simple_test_trial_integrand_functor.hpp"
#include "bempp/fiber/single_component_test_trial_integrand_functor.hpp"
#include "bempp/fiber/surface_div_3d_functor.hpp"
#include "bempp/operators/laplace_operators.hpp"

namespace Bempp {

inline shared_ptr<const ElementaryLocalOperator<double, double>>
identity_operator(const ParameterList &parameterList,
                  const shared_ptr<const Space<double>> &domain,
                  const shared_ptr<const Space<double>> &range,
                  const shared_ptr<const Space<double>> &dualToRange,
                  const std::string &label, int symmetry) {
  return static_pointer_cast<const ElementaryLocalOperator<double, double>>(
      identityOperator<double, double>(parameterList, domain, range,
                                       dualToRange, label, symmetry)
          .abstractOperator());
}

inline shared_ptr<const ElementaryLocalOperator<double, double>>
maxwell_identity_operator(const ParameterList &parameterList,
                          const shared_ptr<const Space<double>> &domain,
                          const shared_ptr<const Space<double>> &range,
                          const shared_ptr<const Space<double>> &dualToRange,
                          const std::string &label, int symmetry) {
  return static_pointer_cast<const ElementaryLocalOperator<double, double>>(
      maxwell3dIdentityOperator<double, double>(parameterList, domain, range,
                                                dualToRange, label, symmetry)
          .abstractOperator());
}

inline shared_ptr<const ElementaryLocalOperator<double, double>>
laplace_beltrami_operator(const ParameterList &parameterList,
                          const shared_ptr<const Space<double>> &domain,
                          const shared_ptr<const Space<double>> &range,
                          const shared_ptr<const Space<double>> &dualToRange,
                          const std::string &label, int symmetry) {
  return static_pointer_cast<const ElementaryLocalOperator<double, double>>(
      laplaceBeltrami3dOperator<double, double>(parameterList, domain, range,
                                                dualToRange, label, symmetry)
          .abstractOperator());
}



// Support operators

inline shared_ptr<const ElementaryLocalOperator<double, double>>
curl_value_local_operator(const shared_ptr<const Space<double>> &domain,
                          const shared_ptr<const Space<double>> &range,
                          const shared_ptr<const Space<double>> &dualToRange,
                          int component) {

  typedef Fiber::ScalarFunctionValueFunctor<double> ValueFunctor;
  typedef Fiber::SurfaceCurl3dFunctor<double> CurlFunctor;
  typedef Fiber::SingleComponentTestTrialIntegrandFunctor<double, double>
      IntegrandFunctor;

  typedef GeneralElementaryLocalOperator<double, double> LocalOp;

  return shared_ptr<const ElementaryLocalOperator<double, double>>(
      new LocalOp(domain, range, dualToRange, "", NO_SYMMETRY, CurlFunctor(),
                  ValueFunctor(), IntegrandFunctor(component, 0)));
}

inline shared_ptr<const ElementaryLocalOperator<double, double>>
value_times_normal_local_operator(const shared_ptr<const Space<double>> &domain,
                          const shared_ptr<const Space<double>> &range,
                          const shared_ptr<const Space<double>> &dualToRange,
                          int component) {

  typedef Fiber::ScalarFunctionValueFunctor<double> ValueFunctor;
  typedef Fiber::ScalarFunctionValueTimesNormalFunctor<double>
      ValueTimesNormalFunctor;
  typedef Fiber::SingleComponentTestTrialIntegrandFunctor<double, double>
      IntegrandFunctor;

  typedef GeneralElementaryLocalOperator<double, double> LocalOp;

  return shared_ptr<const ElementaryLocalOperator<double, double>>(
      new LocalOp(domain, range, dualToRange, "", NO_SYMMETRY, ValueTimesNormalFunctor(),
                  ValueFunctor(), IntegrandFunctor(component, 0)));
}

inline shared_ptr<const ElementaryLocalOperator<double, double>>
vector_value_times_scalar_local_operator(const shared_ptr<const Space<double>> &domain,
                          const shared_ptr<const Space<double>> &range,
                          const shared_ptr<const Space<double>> &dualToRange,
                          int component) {

  typedef Fiber::HdivFunctionValueFunctor<double> VectorValueFunctor;
  typedef Fiber::ScalarFunctionValueFunctor<double> ScalarValueFunctor;
  typedef Fiber::SingleComponentTestTrialIntegrandFunctor<double, double>
      IntegrandFunctor;

  typedef GeneralElementaryLocalOperator<double, double> LocalOp;

  return shared_ptr<const ElementaryLocalOperator<double, double>>(
      new LocalOp(domain, range, dualToRange, "", NO_SYMMETRY, VectorValueFunctor(),
                  ScalarValueFunctor(), IntegrandFunctor(component, 0)));
}

inline shared_ptr<const ElementaryLocalOperator<double, double>>
div_times_scalar_local_operator(const shared_ptr<const Space<double>> &domain,
                          const shared_ptr<const Space<double>> &range,
                          const shared_ptr<const Space<double>> &dualToRange) {

  typedef Fiber::SurfaceDiv3dFunctor<double> DivFunctor;
  typedef Fiber::ScalarFunctionValueFunctor<double> ScalarValueFunctor;
  typedef Fiber::SimpleTestTrialIntegrandFunctor<double, double>
      IntegrandFunctor;

  typedef GeneralElementaryLocalOperator<double, double> LocalOp;

  return shared_ptr<const ElementaryLocalOperator<double, double>>(
      new LocalOp(domain, range, dualToRange, "", NO_SYMMETRY, DivFunctor(),
                  ScalarValueFunctor(), IntegrandFunctor()));
}

}

#endif
