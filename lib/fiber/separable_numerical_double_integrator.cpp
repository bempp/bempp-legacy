// File to be removed once unit tests of SeparableNumericalDoubleIntegrator are in place

#include "separable_numerical_double_integrator.hpp"
#include "../grid/geometry_adapter.hpp"

namespace Fiber
{

#ifdef COMPILE_FOR_FLOAT
template class SeparableNumericalDoubleIntegrator<float, Bempp::GeometryAdapter>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class SeparableNumericalDoubleIntegrator<double, Bempp::GeometryAdapter>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class SeparableNumericalDoubleIntegrator<std::complex<float>, Bempp::GeometryAdapter>;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class SeparableNumericalDoubleIntegrator<std::complex<double>, Bempp::GeometryAdapter>;
#endif

} // namespace Fiber

