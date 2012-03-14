// File to be removed once unit tests of NonseparableNumericalTestKernelTrialIntegrator are in place

#include "nonseparable_numerical_test_kernel_trial_integrator.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/geometry.hpp"

namespace Fiber
{

#ifdef COMPILE_FOR_FLOAT
template class NonseparableNumericalTestKernelTrialIntegrator<float, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class NonseparableNumericalTestKernelTrialIntegrator<double, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class NonseparableNumericalTestKernelTrialIntegrator<std::complex<float>, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class NonseparableNumericalTestKernelTrialIntegrator<std::complex<double>, Bempp::GeometryFactory>;
#endif

} // namespace Fiber

