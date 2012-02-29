// File to be removed once unit tests of NonseparableNumericalDoubleIntegrator are in place

#include "nonseparable_numerical_double_integrator.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/geometry.hpp"

namespace Fiber
{

#ifdef COMPILE_FOR_FLOAT
template class NonseparableNumericalDoubleIntegrator<float, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class NonseparableNumericalDoubleIntegrator<double, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class NonseparableNumericalDoubleIntegrator<std::complex<float>, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class NonseparableNumericalDoubleIntegrator<std::complex<double>, Bempp::GeometryFactory>;
#endif

} // namespace Fiber

