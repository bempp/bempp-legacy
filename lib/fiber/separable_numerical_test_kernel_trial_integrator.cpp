// File to be removed once unit tests of SeparableNumericalTestKernelTrialIntegrator are in place

#include "separable_numerical_test_kernel_trial_integrator.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/geometry.hpp"

namespace Fiber
{

#ifdef COMPILE_FOR_FLOAT
template class SeparableNumericalTestKernelTrialIntegrator<float, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class SeparableNumericalTestKernelTrialIntegrator<double, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class SeparableNumericalTestKernelTrialIntegrator<std::complex<float>, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class SeparableNumericalTestKernelTrialIntegrator<std::complex<double>, Bempp::GeometryFactory>;
#endif

} // namespace Fiber

