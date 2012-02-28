// File to be removed once unit tests of SeparableNumericalDoubleIntegrator are in place

#include "separable_numerical_double_integrator.hpp"
#include "../grid/geometry_factory.hpp"
#include "../grid/geometry.hpp"

namespace Fiber
{

#ifdef COMPILE_FOR_FLOAT
template class SeparableNumericalDoubleIntegrator<float, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class SeparableNumericalDoubleIntegrator<double, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class SeparableNumericalDoubleIntegrator<std::complex<float>, Bempp::GeometryFactory>;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class SeparableNumericalDoubleIntegrator<std::complex<double>, Bempp::GeometryFactory>;
#endif

} // namespace Fiber

