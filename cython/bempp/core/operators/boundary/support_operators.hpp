#ifndef support_operators_hpp
#define support_operators_hpp

#include "bempp/common/boost_make_shared_fwd.hpp"
#include "bempp/common/shared_ptr.hpp"
#include "bempp/assembly/general_elementary_local_operator_imp.hpp"
#include "bempp/assembly/abstract_boundary_operator.hpp"
#include "bempp/assembly/elementary_integral_operator.hpp"
#include "bempp/assembly/elementary_local_operator.hpp"
#include "bempp/assembly/symmetry.hpp"
#include "bempp/fiber/hdiv_function_value_functor.hpp"
#include "bempp/fiber/surface_curl_3d_functor.hpp"
#include "bempp/fiber/scalar_function_value_times_normal_functor.hpp"
#include "bempp/fiber/scalar_function_value_functor.hpp"
#include "bempp/fiber/simple_test_trial_integrand_functor.hpp"
#include "bempp/fiber/single_component_test_trial_integrand_functor.hpp"
#include "bempp/fiber/surface_div_3d_functor.hpp"
#include "bempp/fiber/surface_grad_3d_functor.hpp"
#include "bempp/fiber/hcurl_function_value_functor.hpp"
#include "bempp/operators/laplace_operators.hpp"

namespace Bempp {

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

inline shared_ptr<const ElementaryLocalOperator<double, double>>
div_times_div_local_operator(const shared_ptr<const Space<double>> &domain,
                          const shared_ptr<const Space<double>> &range,
                          const shared_ptr<const Space<double>> &dualToRange) {

  typedef Fiber::SurfaceDiv3dFunctor<double> DivFunctor;
  typedef Fiber::SimpleTestTrialIntegrandFunctor<double, double>
      IntegrandFunctor;

  typedef GeneralElementaryLocalOperator<double, double> LocalOp;

  return shared_ptr<const ElementaryLocalOperator<double, double>>(
      new LocalOp(domain, range, dualToRange, "", NO_SYMMETRY, DivFunctor(),
                  DivFunctor(), IntegrandFunctor()));
}

inline shared_ptr<const ElementaryLocalOperator<double, double>>
grad_times_hcurl_value_local_operator(const shared_ptr<const Space<double>> &domain,
                          const shared_ptr<const Space<double>> &range,
                          const shared_ptr<const Space<double>> &dualToRange) {

  typedef Fiber::SurfaceGrad3dFunctor<double> GradFunctor;
  typedef Fiber::HcurlFunctionValueFunctor<double> HcurlValueFunctor;
  typedef Fiber::SimpleTestTrialIntegrandFunctor<double, double>
      IntegrandFunctor;

  typedef GeneralElementaryLocalOperator<double, double> LocalOp;

  return shared_ptr<const ElementaryLocalOperator<double, double>>(
      new LocalOp(domain, range, dualToRange, "", NO_SYMMETRY, GradFunctor(),
                  HcurlValueFunctor(), IntegrandFunctor()));
}

inline shared_ptr<const ElementaryLocalOperator<double, double>>
hcurl_times_hcurl_value_local_operator(const shared_ptr<const Space<double>> &domain,
                          const shared_ptr<const Space<double>> &range,
                          const shared_ptr<const Space<double>> &dualToRange) {

  typedef Fiber::HcurlFunctionValueFunctor<double> HcurlValueFunctor;
  typedef Fiber::SimpleTestTrialIntegrandFunctor<double, double>
      IntegrandFunctor;

  typedef GeneralElementaryLocalOperator<double, double> LocalOp;

  return shared_ptr<const ElementaryLocalOperator<double, double>>(
      new LocalOp(domain, range, dualToRange, "", NO_SYMMETRY, HcurlValueFunctor(),
                  HcurlValueFunctor(), IntegrandFunctor()));
}



}

#endif
