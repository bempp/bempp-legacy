// File to be removed once unit tests are in place

#include "standard_integration_manager_factory_2d.hpp"
#include "../grid/geometry_adapter.hpp"

namespace Fiber
{

#ifdef COMPILE_FOR_FLOAT
template class StandardIntegrationManagerFactory2D<float, Bempp::GeometryAdapter>;
#endif
#ifdef COMPILE_FOR_DOUBLE
template class StandardIntegrationManagerFactory2D<double, Bempp::GeometryAdapter>;
#endif
#ifdef COMPILE_FOR_COMPLEX_FLOAT
#include <complex>
template class StandardIntegrationManagerFactory2D<std::complex<float>, Bempp::GeometryAdapter>;
#endif
#ifdef COMPILE_FOR_COMPLEX_DOUBLE
#include <complex>
template class StandardIntegrationManagerFactory2D<std::complex<double>, Bempp::GeometryAdapter>;
#endif

} // namespace Fiber
